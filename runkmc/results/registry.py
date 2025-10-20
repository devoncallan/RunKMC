from contextlib import contextmanager
from pathlib import Path
from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional
import sqlite3
import hashlib
import json

import pandas as pd


@dataclass
class SimulationRecord:

    id: int
    timestamp: str
    sim_id: str
    run_index: int
    completed: int
    input_hash: str
    _base_dir: Path

    @staticmethod
    def sql_fields_list(no_id: bool = False) -> List[str]:
        fields = ["id", "timestamp", "sim_id", "run_index", "completed", "input_hash"]
        return fields[1:] if no_id else fields

    @staticmethod
    def sql_fields_str(no_id: bool = False) -> str:
        return ", ".join(SimulationRecord.sql_fields_list(no_id=no_id))

    @property
    def dir(self) -> Path:
        return self._base_dir / self.sim_id


class SimulationRegistry:

    def __init__(self, base_dir: str | Path):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

        self.db_path = self.base_dir / "registry.db"
        self.csv_path = self.base_dir / "registry.csv"
        self._init_db()

    def _init_db(self):
        with self.write_op() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS simulations (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    timestamp TEXT NOT NULL,
                    sim_id TEXT NOT NULL,
                    run_index INTEGER NOT NULL,
                    completed INTEGER NOT NULL DEFAULT 0,
                    input_hash TEXT NOT NULL
                )
                """
            )
            conn.execute(
                "CREATE UNIQUE INDEX IF NOT EXISTS unique_input_run ON simulations(input_hash, run_index)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_input_hash ON simulations(input_hash)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_completed ON simulations(completed)"
            )

    @contextmanager
    def write_op(self):
        conn = sqlite3.connect(self.db_path)
        try:
            yield conn
            conn.commit()
        finally:
            conn.close()
            with sqlite3.connect(self.db_path) as conn:
                df = pd.read_sql_query(
                    f"""
                    SELECT {SimulationRecord.sql_fields_str()} FROM simulations
                    ORDER BY timestamp DESC, id DESC
                    """,
                    conn,
                )
            df["dir"] = df["sim_id"].apply(lambda sid: str(self.base_dir / sid))
            df.to_csv(self.csv_path, index=False)

    @staticmethod
    def hash_input(input_file: Path, params: object | None = None) -> str:
        content = input_file.read_bytes()
        param_str = json.dumps(params, sort_keys=True) if params is not None else ""
        combined = content + param_str.encode()
        return hashlib.sha256(combined).hexdigest()

    @staticmethod
    def make_sim_id(input_hash: str, run_index: int) -> str:
        hash_short = input_hash[:8]
        return f"sim_{hash_short}-r{run_index:03d}"

    @staticmethod
    def _next_run_index(conn: sqlite3.Connection, input_hash: str) -> int:
        (n,) = conn.execute(
            """
            SELECT COALESCE(MAX(run_index), -1) + 1 FROM simulations
            WHERE input_hash = ?
            """,
            (input_hash,),
        ).fetchone()
        return int(n)

    def new_record(self, input_hash: str) -> SimulationRecord:

        timestamp = datetime.now().isoformat()
        with self.write_op() as conn:
            run_index = SimulationRegistry._next_run_index(conn, input_hash)
            sim_id = SimulationRegistry.make_sim_id(input_hash, run_index)
            conn.execute(
                f"""
                INSERT INTO simulations ({SimulationRecord.sql_fields_str(no_id=True)})
                VALUES (?, ?, ?, ?, ?)
                """,
                (timestamp, sim_id, run_index, 0, input_hash),
            )
            rowid = conn.execute("SELECT last_insert_rowid()").fetchone()[0]
        return SimulationRecord(
            rowid, timestamp, sim_id, run_index, 0, input_hash, self.base_dir
        )

    def update_completion(
        self, record: SimulationRecord, completed: bool = True
    ) -> None:
        with self.write_op() as conn:
            conn.execute(
                """
                UPDATE simulations
                SET completed = ?
                WHERE id = ?
                """,
                (1 if completed else 0, record.id),
            )

    def _delete(self, row_id: int) -> bool:
        with self.write_op() as conn:
            cur = conn.execute(
                """
                DELETE FROM simulations
                WHERE id = ?
                """,
                (row_id,),
            )
            return cur.rowcount > 0

    def delete(self, record: SimulationRecord) -> bool:
        with self.write_op() as conn:
            cur = conn.execute(
                """
                DELETE FROM simulations
                WHERE id = ?
                """,
                (record.id,),
            )
            return cur.rowcount > 0

    def latest_completed(self, input_hash: str) -> Optional[SimulationRecord]:
        with sqlite3.connect(self.db_path) as conn:
            row = conn.execute(
                f"""
                SELECT {SimulationRecord.sql_fields_str()} FROM simulations
                WHERE input_hash = ? AND completed = 1
                ORDER BY timestamp DESC
                LIMIT 1
                """,
                (input_hash,),
            ).fetchone()
        return SimulationRecord(*row, _base_dir=self.base_dir) if row else None

    def find_all(self, input_hash: str) -> List[SimulationRecord]:
        with sqlite3.connect(self.db_path) as conn:
            rows = conn.execute(
                f"""
                SELECT {SimulationRecord.sql_fields_str()} FROM simulations
                WHERE input_hash = ?
                ORDER BY timestamp DESC
                """,
                (input_hash,),
            ).fetchall()
        return [SimulationRecord(*row, _base_dir=self.base_dir) for row in rows]

    def get_all(self, completed_only: bool = False) -> List[SimulationRecord]:
        with sqlite3.connect(self.db_path) as conn:
            rows = conn.execute(
                f"""
                SELECT {SimulationRecord.sql_fields_str()} FROM simulations
                {'WHERE completed = 1' if completed_only else ''}
                ORDER BY timestamp DESC
                """
            ).fetchall()
        return [SimulationRecord(*row, _base_dir=self.base_dir) for row in rows]

    def to_df(self) -> pd.DataFrame:
        return (
            pd.read_csv(self.csv_path)
            if self.csv_path.exists()
            else pd.DataFrame(columns=SimulationRecord.sql_fields_list())
        )
