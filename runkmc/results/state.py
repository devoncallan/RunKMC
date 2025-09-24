from __future__ import annotations
from pathlib import Path
from typing import Optional, Dict, Any, List
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
import pandas as pd
import yaml


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
    total_conv: NDArray[np.float64]
    unit_counts: Dict[str, NDArray[np.uint64]]
    polymer_counts: Dict[str, NDArray[np.uint64]]

    # Analysis State
    nAvgCL: NDArray[np.float64]
    wAvgCL: NDArray[np.float64]
    dispCL: NDArray[np.float64]
    nAvgMW: NDArray[np.float64]
    wAvgMW: NDArray[np.float64]
    dispMW: NDArray[np.float64]

    nAvgSL: Dict[str, NDArray[np.float64]]
    wAvgSL: Dict[str, NDArray[np.float64]]
    dispSL: Dict[str, NDArray[np.float64]]

    _raw_data: pd.DataFrame

    @staticmethod
    def from_csv(filepath: Path | str, metadata: Metadata) -> StateData:

        df = pd.read_csv(filepath)

        unit_names = metadata.get_unit_names()
        monomer_names = metadata.get_monomer_names()
        polymer_names = metadata.get_polymer_names()

        return StateData(
            iteration=df["Iteration"].to_numpy(np.uint64),
            kmc_step=df["KMC Step"].to_numpy(np.uint64),
            kmc_time=df["KMC Time"].to_numpy(np.float64),
            sim_time=df["Simulation Time"].to_numpy(np.float64),
            sim_time_per_1e6_steps=df["Simulation Time per 1e6 KMC Steps"].to_numpy(
                np.float64
            ),
            NAV=df["NAV"].to_numpy(np.float64),
            unit_convs={
                name: df[f"Conv_{name}"].to_numpy(np.float64) for name in unit_names
            },
            total_conv=df["Conv_Total"].to_numpy(np.float64),
            unit_counts={
                name: df[f"Count_{name}"].to_numpy(np.uint64) for name in unit_names
            },
            polymer_counts={
                name: df[f"Count_{name}"].to_numpy(np.uint64) for name in polymer_names
            },
            nAvgCL=df["nAvgCL"].to_numpy(np.float64),
            wAvgCL=df["wAvgCL"].to_numpy(np.float64),
            dispCL=df["dispCL"].to_numpy(np.float64),
            nAvgMW=df["nAvgMW"].to_numpy(np.float64),
            wAvgMW=df["wAvgMW"].to_numpy(np.float64),
            dispMW=df["dispMW"].to_numpy(np.float64),
            nAvgSL={
                name: df[f"nAvgSL_{name}"].to_numpy(np.float64)
                for name in monomer_names
            },
            wAvgSL={
                name: df[f"wAvgSL_{name}"].to_numpy(np.float64)
                for name in monomer_names
            },
            dispSL={
                name: df[f"dispSL_{name}"].to_numpy(np.float64)
                for name in monomer_names
            },
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
                if col.startswith("MonCount_"):
                    monomer_names.append(col.replace("MonCount_", ""))

        monomer_names = list(set(monomer_names))

        return SequenceData(
            iteration=df["Iteration"].to_numpy(np.uint64),
            kmc_time=df["KMC Time"].to_numpy(np.float64),
            bucket=df["Bucket"].to_numpy(np.uint64),
            monomer_count={
                name: df[f"monCount_{name}"].to_numpy(np.uint64)
                for name in monomer_names
            },
            sequence_count={
                name: df[f"seqCount_{name}"].to_numpy(np.uint64)
                for name in monomer_names
            },
            sequence_length2={
                name: df[f"seqLengths2_{name}"].to_numpy(np.float64)
                for name in monomer_names
            },
            _raw_data=df,
            _monomer_names=monomer_names,
        )

    @staticmethod
    def from_csv(filepath: Path | str, metadata: Metadata) -> SequenceData:

        df = pd.read_csv(filepath)

        monomer_names = metadata.get_monomer_names()

        return SequenceData._from_df(df, monomer_names)

    def get_buckets(self) -> List[int]:
        return sorted(self._raw_data["Bucket"].unique().tolist())

    def get_by_bucket(self, bucket: int) -> SequenceData:

        valid_buckets = self.get_buckets()
        if bucket not in valid_buckets:
            raise ValueError(
                f"Bucket {bucket} not found in sequence data ({min(valid_buckets)}-{max(valid_buckets)})."
            )

        df = self._raw_data[self._raw_data["Bucket"] == bucket]
        df = df.reset_index(drop=True)
        df = df.sort_values("KMC Time")

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

        required_keys = ["run_info", "parameters", "species", "reactions"]
        missing_keys = [key for key in required_keys if key not in data.keys()]
        if missing_keys:
            raise ValueError(
                f"Metadata file {metadata_path} is missing required keys: {missing_keys}"
            )

        return Metadata(
            run_info=data["run_info"],
            species=data["species"],
            reactions=data["reactions"],
            parameters=data["parameters"],
            _metadata_path=metadata_path,
            _raw_data=data,
        )

    def to_dict(self) -> Dict[str, Any]:
        return self._raw_data if self._raw_data is not None else {}

    def get_monomer_names(self) -> List[str]:

        monomer_names = []
        units: List[Dict[str, Any]] = self.species.get("units", [])
        for unit in units:
            if unit["type"] == "M":
                monomer_names.append(unit["name"])

        if len(monomer_names) == 0:
            raise ValueError(f"Metadata file does not contain any monomer information.")

        return monomer_names

    def get_unit_names(self) -> List[str]:

        unit_names = []

        units: List[Dict[str, Any]] = self.species.get("units", [])
        if len(units) == 0:
            raise ValueError(f"Metadata file does not contain any unit information.")

        for unit in units:
            unit_names.append(unit["name"])

        return unit_names

    def get_polymer_names(self) -> List[str]:

        polymer_names = []

        polymers: List[Dict[str, Any]] = self.species.get("polymers", [])
        if len(polymers) == 0:
            raise ValueError(f"Metadata file does not contain any polymer information.")

        for polymer in polymers:
            polymer_names.append(polymer["name"])

        return polymer_names
