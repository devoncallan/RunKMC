from __future__ import annotations
from pathlib import Path
from typing import Optional, Dict, Any, List
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
import pandas as pd
import yaml

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

        df = pd.read_csv(filepath)

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

        monomer_names = list(set(monomer_names))

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
            df = pd.read_csv(filepath)
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
        nAvgSL = {}
        for name in self._monomer_names:
            nAvgSL[name] = np.divide(
                self.monomer_count[name],
                self.sequence_count[name],
                where=self.sequence_count[name] != 0,
            )
        return nAvgSL

    @property
    def wAvgSL(self) -> Dict[str, NDArray[np.float64]]:
        wAvgSL = {}
        for name in self._monomer_names:
            wAvgSL[name] = np.divide(
                self.sequence_length2[name],
                self.sequence_count[name],
                where=self.sequence_count[name] != 0,
            )
        return wAvgSL
