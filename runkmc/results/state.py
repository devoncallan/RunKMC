from __future__ import annotations
from pathlib import Path
from typing import Optional, Dict, Any, List
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
import pandas as pd
import yaml

from ..core import C


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
    def from_csv(filepath: Path | str, metadata: Metadata) -> StateData:

        df = pd.read_csv(filepath)

        unit_names = metadata.get_unit_names()
        monomer_names = metadata.get_monomer_names()
        polymer_names = metadata.get_polymer_names()

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
    def from_csv(filepath: Path | str, metadata: Metadata) -> SequenceData:

        try:
            df = pd.read_csv(filepath)
            monomer_names = metadata.get_monomer_names()
            return SequenceData._from_df(df, monomer_names)
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


@dataclass
class Metadata:
    run_info: Dict[str, Any]
    species: Dict[str, Any]
    reactions: Dict[str, Any]
    parameters: Dict[str, Any]

    _metadata_path: Optional[Path] = None
    _raw_data: Optional[Dict[str, Any]] = None

    @staticmethod
    def load(metadata_path: Path | str) -> Metadata:

        metadata_path = Path(metadata_path)
        with open(metadata_path, "r") as file:
            data = yaml.safe_load(file)

        assert (
            data is not None
        ), f"Metadata file {metadata_path} is empty or invalid YAML."
        assert isinstance(
            data, dict
        ), f"Metadata file {metadata_path} does not contain a valid YAML dictionary."

        required_keys = [
            C.io.RUN_INFO_SECTION,
            C.io.PARAMETERS_SECTION,
            C.io.SPECIES_SECTION,
            C.io.REACTIONS_SECTION,
        ]
        missing_keys = [key for key in required_keys if key not in data.keys()]
        if missing_keys:
            raise ValueError(
                f"Metadata file {metadata_path} is missing required keys: {missing_keys}"
            )

        return Metadata(
            run_info=data[C.io.RUN_INFO_SECTION],
            species=data[C.io.SPECIES_SECTION],
            reactions=data[C.io.REACTIONS_SECTION],
            parameters=data[C.io.PARAMETERS_SECTION],
            _metadata_path=metadata_path,
            _raw_data=data,
        )

    def to_dict(self) -> Dict[str, Any]:
        return self._raw_data if self._raw_data is not None else {}

    def get_monomer_names(self) -> List[str]:

        monomer_names = []
        units: List[Dict[str, Any]] = self.species.get(C.io.UNITS_KEY, [])
        for unit in units:
            if unit[C.io.TYPE_KEY] == "M":
                monomer_names.append(unit[C.io.NAME_KEY])

        if len(monomer_names) == 0:
            raise ValueError(f"Metadata file does not contain any monomer information.")

        return monomer_names

    def get_unit_names(self) -> List[str]:

        unit_names = []

        units: List[Dict[str, Any]] = self.species.get(C.io.UNITS_KEY, [])
        if len(units) == 0:
            raise ValueError(f"Metadata file does not contain any unit information.")

        for unit in units:
            unit_names.append(unit[C.io.NAME_KEY])

        return unit_names

    def get_polymer_names(self) -> List[str]:

        polymer_names = []

        polymers: List[Dict[str, Any]] = self.species.get(C.io.POLYMERS_KEY, [])
        if len(polymers) == 0:
            raise ValueError(f"Metadata file does not contain any polymer information.")

        for polymer in polymers:
            polymer_names.append(polymer[C.io.NAME_KEY])

        return polymer_names
