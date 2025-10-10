from pathlib import Path
from dataclasses import dataclass
import sqlite3


@dataclass
class SimulationRecord:
    sim_id: str
    timestamp: str


class SimulationRegistry:

    def __init__(self, output_dir: str | Path):
        self.output_dir = Path(output_dir)
        self.db_path = self.output_dir / "registry.db"
        self._init_db()

    def _init_db(self):
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS simulations (
                    hash TEXT PRIMARY KEY,
                    output_dir TEXT,
                    timestamp TEXT,
                    completed INTEGER
                )
                """
            )
            conn.execute("CREATE INDEX IF NOT EXISTS idx_hash ON simulations(hash)")

    def find_by_hash(self, input_hash: str) -> SimulationRecord | None:
        with sqlite3.connect(self.db_path) as conn:
            row = conn.execute(
                "SELECT * FROM simulations WHERE hash = ?", (input_hash,)
            ).fetchone()
            return SimulationRecord(*row) if row else None
