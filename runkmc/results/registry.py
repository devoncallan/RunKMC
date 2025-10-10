from contextlib import contextmanager
from pathlib import Path
from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional
import sqlite3
import random
import hashlib

import pandas as pd
import friendlywords as fw


@dataclass
class SimulationRecord:
    """Record of a simulation in the registry database.

    Attributes:
        input_hash: SHA256 hash of the parsed input file
        sim_id: Unique identifier for the simulation
        timestamp: ISO format timestamp of when simulation was registered
        completed: Whether the simulation completed successfully (1) or not (0)
    """

    id: int
    timestamp: str

    sim_id: str
    run_index: int
    completed: int
    input_hash: str

    _base_dir: Path  # not included in db

    @staticmethod
    def sql_fields_list(no_id: bool = False) -> List[str]:
        fields = ["id", "timestamp", "sim_id", "run_index", "completed", "input_hash"]
        return fields[1:] if no_id else fields

    @staticmethod
    def sql_fields_str(no_id: bool = False) -> str:
        return ", ".join(SimulationRecord.sql_fields_list(no_id=no_id))

    def to_sql_row(self, no_id: bool = False) -> tuple:
        return (
            self.id,
            self.timestamp,
            self.sim_id,
            self.run_index,
            self.completed,
            self.input_hash,
        )

    @property
    def dir(self) -> Path:
        return self._base_dir / self.sim_id


class SimulationRegistry:
    """Registry for tracking simulations and enabling result caching.

    Uses SQLite database to store mappings between input file hashes
    and their corresponding output directories. This enables:
    - Automatic detection of duplicate simulation requests
    - Fast lookup of existing results
    - Tracking of simulation completion status
    """

    def __init__(self, base_dir: str | Path):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

        self.db_path = self.base_dir / "registry.db"
        self.csv_path = self.base_dir / "registry.csv"
        self._init_db()

    def _init_db(self):
        """Initialize the registry database with required schema."""
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
        """Context manager for writing operations with automatic CSV export."""
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

    # ------------ Helpers ------------
    @staticmethod
    def hash_file(filepath: Path) -> str:
        """Compute SHA256 hash of a file."""
        return hashlib.sha256(filepath.read_bytes()).hexdigest()

    @staticmethod
    def generate_base_name(hash_str: str) -> str:
        # Save and restore random state to avoid global side effects
        rs = random.getstate()
        random.seed(hash_str)
        base = str(fw.generate("po", separator="_"))
        random.setstate(rs)
        return base

    def _make_sim_id(self, input_hash: str, run_index: int) -> str:
        base_name = self.generate_base_name(input_hash)
        return f"{base_name}_{input_hash[:4]}-r{run_index:03d}"

    def _next_run_index(self, conn: sqlite3.Connection, input_hash: str) -> int:
        (n,) = conn.execute(
            """
            SELECT COALESCE(MAX(run_index), -1) + 1 FROM simulations
            WHERE input_hash = ?
            """,
            (input_hash,),
        ).fetchone()
        return int(n)

    # ------------ API ------------

    def insert(self, input_hash: str) -> SimulationRecord:
        """Register a new simulation in the database.
        Args:
            input_hash: SHA256 hash of the parsed input file
        Returns:
            SimulationRecord of the newly inserted simulation
        """

        timestamp = datetime.now().isoformat()
        with self.write_op() as conn:
            run_index = self._next_run_index(conn, input_hash)
            sim_id = self._make_sim_id(input_hash, run_index)
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

    def update_completion(self, row_id: int, completed: bool = True) -> None:
        """Mark a simulation as completed or not completed.
        Args:
            row_id: Row ID of the simulation record to update
            completed: Completion status (default: True)
        """
        with self.write_op() as conn:
            conn.execute(
                """
                UPDATE simulations
                SET completed = ? 
                WHERE id = ?
                """,
                (1 if completed else 0, row_id),
            )

    def delete(self, row_id: int) -> bool:
        """Remove a simulation from the registry.
        Args:
            row_id: Row ID of the simulation record to delete
        Returns:
            True if a record was deleted, False if not found
        """
        with self.write_op() as conn:
            cur = conn.execute(
                """
                DELETE FROM simulations
                WHERE id = ?
                """,
                (row_id,),
            )
            return cur.rowcount > 0

    def latest_completed(self, input_hash: str) -> Optional[SimulationRecord]:
        """Get the latest completed simulation for a given input hash.
        Args:
            input_hash: SHA256 hash of the parsed input file
        Returns:
            SimulationRecord if found, None otherwise
        """

        with sqlite3.connect(self.db_path) as conn:
            row = conn.execute(
                f"""
                SELECT {SimulationRecord.sql_fields_str()} FROM simulations
                WHERE input_hash = ? AND completed = 1
                ORDER BY run_index DESC
                LIMIT 1
                """,
                (input_hash,),
            ).fetchone()
        return SimulationRecord(*row, self.base_dir) if row else None

    def find_all(self, input_hash: str) -> List[SimulationRecord]:
        """Find all simulations for a given input hash.
        Args:
            input_hash: SHA256 hash of the parsed input file
        Returns:
            List of SimulationRecord objects
        """
        with sqlite3.connect(self.db_path) as conn:
            rows = conn.execute(
                f"""
                SELECT {SimulationRecord.sql_fields_str()} FROM simulations
                WHERE input_hash = ?
                ORDER BY run_index DESC
                """,
                (input_hash,),
            ).fetchall()
        return [SimulationRecord(*row, self.base_dir) for row in rows]

    def get_all(self, completed_only: bool = False) -> List[SimulationRecord]:
        """Get all registered simulations.
        Args:
            completed_only: If True, only return completed simulations
        Returns:
            List of SimulationRecord objects
        """
        with sqlite3.connect(self.db_path) as conn:
            rows = conn.execute(
                f"""
                SELECT {SimulationRecord.sql_fields_str()} FROM simulations
                {'WHERE completed = 1' if completed_only else ''}
                ORDER BY timestamp DESC
                """
            ).fetchall()
        return [SimulationRecord(*row, self.base_dir) for row in rows]

    def to_df(self) -> pd.DataFrame:
        """Export the registry to a pandas DataFrame.
        Returns:
            DataFrame containing all simulation records
        """

        return (
            pd.read_csv(self.csv_path)
            if self.csv_path.exists()
            else pd.DataFrame(columns=SimulationRecord.sql_fields_list())
        )
