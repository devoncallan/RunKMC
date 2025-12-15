from __future__ import annotations
from pathlib import Path
from typing import Optional, Dict, List, Any, Tuple
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
import pandas as pd

from ..core import C
from ..core.species import SpeciesRegistry


@dataclass
class StateData:

    # KMC State
    iteration: NDArray[np.uint64]
    kmc_step: NDArray[np.uint64]
    kmc_time: NDArray[np.float64]
    sim_time: NDArray[np.float64]
    sim_time_per_1e6_steps: NDArray[np.float64]
    NAV: NDArray[np.float64]

    # Species State
    unit_convs: Dict[str, NDArray[np.float64]]
    monomer_conv: NDArray[np.float64]
    unit_counts: Dict[str, NDArray[np.uint64]]
    polymer_counts: Dict[str, NDArray[np.uint64]]

    # Analysis State
    nAvgCL: NDArray[np.float64]
    wAvgCL: NDArray[np.float64]
    dispCL: NDArray[np.float64]
    nAvgMW: NDArray[np.float64]
    wAvgMW: NDArray[np.float64]
    dispMW: NDArray[np.float64]

    nAvgSL: Dict[str, NDArray[np.float64]] | None
    wAvgSL: Dict[str, NDArray[np.float64]] | None
    dispSL: Dict[str, NDArray[np.float64]] | None

    _raw_data: pd.DataFrame

    @staticmethod
    def from_csv(filepath: Path | str, species: SpeciesRegistry) -> StateData:

        df = pd.read_csv(filepath)  # type: ignore

        unit_names = species.get_unit_names()
        monomer_names = species.get_monomer_names()
        polymer_names = species.get_polymer_names()

        no_sequence = len(monomer_names) <= 1

        return StateData(
            iteration=df[C.state.ITERATION_KEY].to_numpy(np.uint64),
            kmc_step=df[C.state.KMC_STEP_KEY].to_numpy(np.uint64),
            kmc_time=df[C.state.KMC_TIME_KEY].to_numpy(np.float64),
            sim_time=df[C.state.SIM_TIME_KEY].to_numpy(np.float64),
            sim_time_per_1e6_steps=df[C.state.SIM_TIME_PER_1E6_STEPS_KEY].to_numpy(
                np.float64
            ),
            NAV=df[C.state.NAV_KEY].to_numpy(np.float64),
            unit_convs={
                name: df[C.state.CONV_PREFIX + name].to_numpy(np.float64)
                for name in unit_names
            },
            monomer_conv=df[C.state.CONV_PREFIX + C.state.MONOMER].to_numpy(np.float64),
            unit_counts={
                name: df[C.state.COUNT_PREFIX + name].to_numpy(np.uint64)
                for name in unit_names
            },
            polymer_counts={
                name: df[C.state.COUNT_PREFIX + name].to_numpy(np.uint64)
                for name in polymer_names
            },
            nAvgCL=df[C.state.NAVGCL_KEY].to_numpy(np.float64),
            wAvgCL=df[C.state.WAVGCL_KEY].to_numpy(np.float64),
            dispCL=df[C.state.DISPCL_KEY].to_numpy(np.float64),
            nAvgMW=df[C.state.NAVGMW_KEY].to_numpy(np.float64),
            wAvgMW=df[C.state.WAVGMW_KEY].to_numpy(np.float64),
            dispMW=df[C.state.DISPMW_KEY].to_numpy(np.float64),
            nAvgSL=(
                {
                    name: df[C.state.NAVGSL_PREFIX + name].to_numpy(np.float64)
                    for name in monomer_names
                }
                if not no_sequence
                else None
            ),
            wAvgSL=(
                {
                    name: df[C.state.WAVGSL_PREFIX + name].to_numpy(np.float64)
                    for name in monomer_names
                }
                if not no_sequence
                else None
            ),
            dispSL=(
                {
                    name: df[C.state.DISP_SL_PREFIX + name].to_numpy(np.float64)
                    for name in monomer_names
                }
                if not no_sequence
                else None
            ),
            _raw_data=df,
        )


@dataclass
class SequenceData:

    iteration: NDArray[np.uint64]
    kmc_time: NDArray[np.float64]
    bucket: NDArray[np.uint64]
    monomer_count: Dict[str, NDArray[np.uint64]]
    sequence_count: Dict[str, NDArray[np.uint64]]
    sequence_length2: Dict[str, NDArray[np.float64]]

    _raw_data: pd.DataFrame
    _monomer_names: List[str]

    @staticmethod
    def _from_df(df: pd.DataFrame, monomer_names: Optional[List[str]]) -> SequenceData:

        if monomer_names is None or len(monomer_names) == 0:
            monomer_names = []
            for col in df.columns:
                if col.startswith(C.state.MONCOUNT_PREFIX):
                    monomer_names.append(col.replace(C.state.MONCOUNT_PREFIX, ""))

        monomer_names = list(dict.fromkeys(monomer_names))

        return SequenceData(
            iteration=df[C.state.ITERATION_KEY].to_numpy(np.uint64),
            kmc_time=df[C.state.KMC_TIME_KEY].to_numpy(np.float64),
            bucket=df[C.state.BUCKET_KEY].to_numpy(np.uint64),
            monomer_count={
                name: df[C.state.MONCOUNT_PREFIX + name].to_numpy(np.uint64)
                for name in monomer_names
            },
            sequence_count={
                name: df[C.state.SEQCOUNT_PREFIX + name].to_numpy(np.uint64)
                for name in monomer_names
            },
            sequence_length2={
                name: df[C.state.SEQLEN2_PREFIX + name].to_numpy(np.float64)
                for name in monomer_names
            },
            _raw_data=df,
            _monomer_names=monomer_names,
        )

    @staticmethod
    def from_csv(filepath: Path | str, species: SpeciesRegistry) -> SequenceData:

        try:
            df = pd.read_csv(filepath)  # type: ignore
            return SequenceData._from_df(df, species.get_monomer_names())
        except Exception as e:
            raise ValueError(f"Error loading sequence data from {filepath}: {e}")

    def get_buckets(self) -> List[int]:
        return sorted(self._raw_data[C.state.BUCKET_KEY].unique().tolist())

    def get_by_bucket(self, bucket: int) -> SequenceData:

        valid_buckets = self.get_buckets()
        if bucket not in valid_buckets:
            raise ValueError(
                f"Bucket {bucket} not found in sequence data ({min(valid_buckets)}-{max(valid_buckets)})."
            )

        df = self._raw_data[self._raw_data[C.state.BUCKET_KEY] == bucket]
        df = df.reset_index(drop=True)
        df = df.sort_values(C.state.KMC_TIME_KEY)

        return SequenceData._from_df(df, self._monomer_names)

    @property
    def nAvgSL(self) -> Dict[str, NDArray[np.float64]]:
        nAvgSL: Dict[str, NDArray[np.float64]] = {}
        for name in self._monomer_names:
            nAvgSL[name] = np.divide(
                self.monomer_count[name],
                self.sequence_count[name],
                where=self.sequence_count[name] != 0,
            )
        return nAvgSL

    @property
    def wAvgSL(self) -> Dict[str, NDArray[np.float64]]:
        wAvgSL: Dict[str, NDArray[np.float64]] = {}
        for name in self._monomer_names:
            wAvgSL[name] = np.divide(
                self.sequence_length2[name],
                self.sequence_count[name],
                where=self.sequence_count[name] != 0,
            )
        return wAvgSL


@dataclass
class ChainHistogramData:

    @dataclass
    class ChainHistogram:
        iteration: int
        kmc_time: float
        data: pd.DataFrame
        monomer_names: List[str]
        has_sequence_stats: bool

        @property
        def chain_count(self) -> pd.Series:
            return self.data[C.state.CHAINCOUNT_KEY]

        def bin_counts(self, monomer: str) -> pd.Series:
            column = f"{C.state.BIN_MONCOUNT_PREFIX}{monomer}"
            if column not in self.data.columns:
                raise KeyError(f"Bin counts for monomer '{monomer}' not found.")
            return self.data[column]

        def total_sequence_counts(self, monomer: str) -> Optional[pd.Series]:
            if not self.has_sequence_stats:
                return None
            column = f"{C.state.TOTAL_SEQCOUNT_PREFIX}{monomer}"
            return self.data[column] if column in self.data.columns else None

        def total_sequence_lengths2(self, monomer: str) -> Optional[pd.Series]:
            if not self.has_sequence_stats:
                return None
            column = f"{C.state.TOTAL_SEQLEN2_PREFIX}{monomer}"
            return self.data[column] if column in self.data.columns else None

    blocks: List[ChainHistogram]
    monomer_names: List[str]
    has_sequence_stats: bool

    @property
    def iterations(self) -> NDArray[np.uint64]:
        return np.asarray([block.iteration for block in self.blocks], dtype=np.uint64)

    @property
    def times(self) -> NDArray[np.float64]:
        return np.asarray([block.kmc_time for block in self.blocks], dtype=np.float64)

    def latest(self) -> ChainHistogram:
        if not self.blocks:
            raise ValueError("No histogram blocks available.")
        return self.blocks[-1]

    def cumulative_until(self, index: int) -> ChainHistogram:
        if index < 0:
            index += len(self.blocks)
        if index < 0 or index >= len(self.blocks):
            raise IndexError("Histogram block index out of range.")
        return self.blocks[index]

    def concat(self) -> pd.DataFrame:
        if not self.blocks:
            return pd.DataFrame()
        frames: List[pd.DataFrame] = []
        for block in self.blocks:
            frame = block.data.copy()
            frame[C.state.ITERATION_KEY] = block.iteration
            frame[C.state.KMC_TIME_KEY] = block.kmc_time
            frames.append(frame)
        return pd.concat(frames, ignore_index=True)

    def chain_lengths(self, *, block: int = -1) -> np.ndarray:
        if not self.blocks:
            return np.asarray([], dtype=float)
        block = len(self.blocks) + block if block < 0 else block
        if block < 0 or block >= len(self.blocks):
            raise IndexError("Histogram block index out of range.")

        columns = [
            f"{C.state.BIN_MONCOUNT_PREFIX}{name}"
            for name in self.monomer_names
            if f"{C.state.BIN_MONCOUNT_PREFIX}{name}" in self.blocks[block].data.columns
        ]
        if not columns:
            return np.asarray([], dtype=float)

        lengths = self.blocks[block].data[columns].sum(axis=1).to_numpy(dtype=float)
        return lengths

    def chain_counts(self, *, block: int = -1) -> np.ndarray:
        if not self.blocks:
            return np.asarray([], dtype=float)
        block = len(self.blocks) + block if block < 0 else block
        if block < 0 or block >= len(self.blocks):
            raise IndexError("Histogram block index out of range.")

        counts = self.blocks[block].chain_count.to_numpy(dtype=float)
        return counts

    @staticmethod
    def load(filepath: Path | str, species: SpeciesRegistry) -> Optional["ChainHistogramData"]:
        records = read_histogram_data(filepath)
        if not records:
            return None

        blocks: List[ChainHistogramData.ChainHistogram] = []
        monomer_names: List[str] = []
        has_sequence_stats = False

        for header, df in records:
            iteration = int(header[C.state.ITERATION_KEY])
            kmc_time = float(header[C.state.KMC_TIME_KEY])
            frame = df.reset_index(drop=True)

            if not monomer_names:
                monomer_names = [
                    column.replace(C.state.BIN_MONCOUNT_PREFIX, "")
                    for column in frame.columns
                    if column.startswith(C.state.BIN_MONCOUNT_PREFIX)
                ]
                if not monomer_names:
                    monomer_names = species.get_monomer_names()

                has_sequence_stats = (
                    len(monomer_names) > 1
                    and all(
                        f"{C.state.TOTAL_SEQCOUNT_PREFIX}{name}" in frame.columns
                        and f"{C.state.TOTAL_SEQLEN2_PREFIX}{name}" in frame.columns
                        for name in monomer_names
                    )
                )

            blocks.append(
                ChainHistogramData.ChainHistogram(
                    iteration=iteration,
                    kmc_time=kmc_time,
                    data=frame,
                    monomer_names=monomer_names,
                    has_sequence_stats=has_sequence_stats,
                )
            )

        return ChainHistogramData(
            blocks=blocks,
            monomer_names=monomer_names,
            has_sequence_stats=has_sequence_stats,
        )


def parse_block_header(line: str) -> Dict[str, Any]:

    parts: Dict[str, str] = dict(kv.split("=", 1) for kv in line[1:].split(",") if "=" in kv)
    try:
        iteration = int(parts[C.state.ITERATION_KEY])
        kmc_time = float(parts[C.state.KMC_TIME_KEY])
        bins = int(parts[C.state.BINS_KEY])
    except KeyError as exc:
        raise ValueError(f"Missing histogram header field: {exc}") from exc
    return {
        C.state.ITERATION_KEY: iteration,
        C.state.KMC_TIME_KEY: kmc_time,
        C.state.BINS_KEY: bins,
    }


def read_histogram_data(
    path: Path | str,
) -> List[Tuple[Dict[str, Any], pd.DataFrame]]:

    data = []
    expected_columns: List[str] = []
    with Path(path).open("r", encoding="utf-8") as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("#"):
                header = parse_block_header(stripped)

                df = pd.read_csv(fh, nrows=int(header[C.state.BINS_KEY]), header=0)
                current_columns = list(df.columns)
                if expected_columns:
                    if current_columns != expected_columns:
                        raise ValueError("Inconsistent chain histogram columns.")
                else:
                    expected_columns = current_columns

                assert len(df) == int(
                    header[C.state.BINS_KEY]
                ), "Histogram bin count mismatch."
                data.append((header, df))

    return data
