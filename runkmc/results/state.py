from __future__ import annotations
import io
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
    terminated_chain_count: NDArray[np.uint64]

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
            terminated_chain_count=df[C.state.TERMINATED_CHAIN_COUNT_KEY].to_numpy(np.uint64),
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
class ChainInterval:
    """Chains terminated during a single output interval."""
    iteration: int
    kmc_time: float
    data: pd.DataFrame
    monomer_names: List[str]
    has_sequence_stats: bool

    def monomer_counts(self, monomer: str) -> NDArray[np.uint64]:
        col = C.state.MONCOUNT_PREFIX + monomer
        return self.data[col].to_numpy(dtype=np.uint64)

    def sequence_counts(self, monomer: str) -> Optional[NDArray[np.uint64]]:
        if not self.has_sequence_stats:
            return None
        col = C.state.SEQCOUNT_PREFIX + monomer
        return self.data[col].to_numpy(dtype=np.uint64) if col in self.data.columns else None

    def sequence_lengths2(self, monomer: str) -> Optional[NDArray[np.float64]]:
        if not self.has_sequence_stats:
            return None
        col = C.state.SEQLEN2_PREFIX + monomer
        return self.data[col].to_numpy(dtype=np.float64) if col in self.data.columns else None

    def chain_lengths(self) -> NDArray[np.uint64]:
        cols = [C.state.MONCOUNT_PREFIX + m for m in self.monomer_names]
        return self.data[cols].sum(axis=1).to_numpy(dtype=np.uint64)


@dataclass
class ChainRecordData:
    """Per-chain records from chain_stats.csv, optionally split by output interval."""

    data: pd.DataFrame
    monomer_names: List[str]
    has_sequence_stats: bool
    intervals: Optional[List[ChainInterval]]

    @property
    def chain_lengths(self) -> NDArray[np.uint64]:
        cols = [C.state.MONCOUNT_PREFIX + m for m in self.monomer_names]
        return self.data[cols].sum(axis=1).to_numpy(dtype=np.uint64)

    def monomer_counts(self, monomer: str) -> NDArray[np.uint64]:
        col = C.state.MONCOUNT_PREFIX + monomer
        return self.data[col].to_numpy(dtype=np.uint64)

    def sequence_counts(self, monomer: str) -> Optional[NDArray[np.uint64]]:
        if not self.has_sequence_stats:
            return None
        col = C.state.SEQCOUNT_PREFIX + monomer
        return self.data[col].to_numpy(dtype=np.uint64) if col in self.data.columns else None

    def sequence_lengths2(self, monomer: str) -> Optional[NDArray[np.float64]]:
        if not self.has_sequence_stats:
            return None
        col = C.state.SEQLEN2_PREFIX + monomer
        return self.data[col].to_numpy(dtype=np.float64) if col in self.data.columns else None

    def get_interval(self, index: int) -> ChainInterval:
        if self.intervals is None:
            raise ValueError("No interval data available (load with state_data to enable).")
        if index < 0:
            index += len(self.intervals)
        if index < 0 or index >= len(self.intervals):
            raise IndexError("Interval index out of range.")
        return self.intervals[index]

    def cumulative_until(self, index: int) -> ChainInterval:
        """Return a ChainInterval aggregating all chains up to and including index."""
        if self.intervals is None:
            raise ValueError("No interval data available (load with state_data to enable).")
        if index < 0:
            index += len(self.intervals)
        if index < 0 or index >= len(self.intervals):
            raise IndexError("Interval index out of range.")
        combined = pd.concat(
            [iv.data for iv in self.intervals[: index + 1]], ignore_index=True
        )
        last = self.intervals[index]
        return ChainInterval(
            iteration=last.iteration,
            kmc_time=last.kmc_time,
            data=combined,
            monomer_names=self.monomer_names,
            has_sequence_stats=self.has_sequence_stats,
        )

    @staticmethod
    def load(
        filepath: Path | str,
        species: SpeciesRegistry,
        state_data: Optional["StateData"] = None,
    ) -> Optional["ChainRecordData"]:
        df = pd.read_csv(filepath)  # type: ignore
        if df.empty:
            return None

        monomer_names = [
            col.replace(C.state.MONCOUNT_PREFIX, "")
            for col in df.columns
            if col.startswith(C.state.MONCOUNT_PREFIX)
        ]
        if not monomer_names:
            monomer_names = species.get_monomer_names()

        has_sequence_stats = len(monomer_names) > 1 and all(
            C.state.SEQCOUNT_PREFIX + m in df.columns for m in monomer_names
        )

        intervals: Optional[List[ChainInterval]] = None
        if state_data is not None:
            intervals = []
            offset = 0
            counts = state_data.terminated_chain_count
            iterations = state_data.iteration
            times = state_data.kmc_time
            for i, n in enumerate(counts):
                n = int(n)
                slice_df = df.iloc[offset : offset + n].reset_index(drop=True)
                intervals.append(
                    ChainInterval(
                        iteration=int(iterations[i]),
                        kmc_time=float(times[i]),
                        data=slice_df,
                        monomer_names=monomer_names,
                        has_sequence_stats=has_sequence_stats,
                    )
                )
                offset += n

        return ChainRecordData(
            data=df,
            monomer_names=monomer_names,
            has_sequence_stats=has_sequence_stats,
            intervals=intervals,
        )


@dataclass
class SegmentHistogramBlock:
    iteration: int
    kmc_time: float
    data: pd.DataFrame
    monomer_names: List[str]

    def histogram(
        self, monomer_id: str
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        column = C.state.SEGMENT_COUNT_PREFIX + monomer_id
        if column not in self.data.columns:
            raise KeyError(f"Segment histogram for monomer '{monomer_id}' not found.")
        lengths = self.data[C.state.SEGMENT_LENGTH_KEY].to_numpy(dtype=np.float64)
        counts = self.data[column].to_numpy(dtype=np.float64)
        return lengths, counts


@dataclass
class SegmentHistogramData:
    blocks: List[SegmentHistogramBlock]
    monomer_names: List[str]

    @property
    def iterations(self) -> NDArray[np.uint64]:
        return np.asarray([block.iteration for block in self.blocks], dtype=np.uint64)

    @property
    def times(self) -> NDArray[np.float64]:
        return np.asarray([block.kmc_time for block in self.blocks], dtype=np.float64)

    def latest(self) -> SegmentHistogramBlock:
        if not self.blocks:
            raise ValueError("No segment histogram blocks available.")
        return self.blocks[-1]

    def get_histogram(
        self, monomer_id: str, block: int = -1
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        if not self.blocks:
            return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)
        if monomer_id not in self.monomer_names:
            raise ValueError(
                f"Monomer '{monomer_id}' not available in segment histogram."
            )
        block_idx = len(self.blocks) + block if block < 0 else block
        if block_idx < 0 or block_idx >= len(self.blocks):
            raise IndexError("Segment histogram block index out of range.")
        return self.blocks[block_idx].histogram(monomer_id)

    @staticmethod
    def load(
        filepath: Path | str, species: SpeciesRegistry
    ) -> Optional["SegmentHistogramData"]:
        records = read_histogram_data(filepath)
        if not records:
            return None

        blocks: List[SegmentHistogramBlock] = []
        monomer_names: List[str] = []

        for header, df in records:
            if df.empty or C.state.SEGMENT_LENGTH_KEY not in df.columns:
                continue
            if not monomer_names:
                monomer_names = [
                    col.replace(C.state.SEGMENT_COUNT_PREFIX, "")
                    for col in df.columns
                    if col.startswith(C.state.SEGMENT_COUNT_PREFIX)
                ]
                if not monomer_names:
                    monomer_names = species.get_monomer_names()

            blocks.append(
                SegmentHistogramBlock(
                    iteration=int(header[C.state.ITERATION_KEY]),
                    kmc_time=float(header[C.state.KMC_TIME_KEY]),
                    data=df.reset_index(drop=True),
                    monomer_names=monomer_names,
                )
            )

        if not blocks:
            return None

        return SegmentHistogramData(blocks=blocks, monomer_names=monomer_names)


def parse_block_header(line: str) -> Dict[str, Any]:

    parts: Dict[str, str] = dict(
        kv.split("=", 1) for kv in line[1:].split(",") if "=" in kv
    )
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
    expected_cols: List[str] | None = None
    with Path(path).open("r", encoding="utf-8") as fh:

        while True:
            line = fh.readline()
            if not line:
                break
            if not line.startswith("#"):
                continue

            header = parse_block_header(line.strip())
            bins = int(header[C.state.BINS_KEY])

            col_line = fh.readline()
            if not col_line:
                raise ValueError("Missing histogram column header line.")

            rows = []
            for _ in range(bins):
                row = fh.readline()
                if not row:
                    raise ValueError(
                        f"EOF while reading Iteration={header[C.state.ITERATION_KEY]} histogram data."
                    )
                rows.append(row)

            df = pd.read_csv(io.StringIO(col_line + "".join(rows)))

            if len(df) != bins:
                raise ValueError(
                    f"Expected {bins} rows for Iteration={header[C.state.ITERATION_KEY]} histogram, got {len(df)}."
                )

            cols = list(df.columns)
            if expected_cols is None:
                expected_cols = cols
            elif cols != expected_cols:
                raise ValueError(
                    f"Inconsistent histogram columns. Expected {expected_cols}, got {cols}."
                )

            data.append((header, df))
    return data
