"""Parameterized positive-sequence algebraic network model for IEEE 9-bus.

The checked-in load-flow workbook snapshot provides branch topology and
operating-point flows, but not branch R/X/charging parameters. This module
therefore requires those physical parameters explicitly; it never infers them
silently from a single power-flow point.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
BUS_IDS = tuple(range(1, 10))
PARAMETER_COLUMNS = ("r_pu", "x_pu", "b_pu", "tap_pu", "shift_deg")


class IEEE9BusNetwork:
    """Positive-sequence bus network with pi branches and constant-Z loads.

    ``branch_parameters`` must contain one row per ``name`` from
    ``IEEE9_connections_loadflow.csv`` and physical columns `r_pu`, `x_pu`,
    `b_pu` (total line charging susceptance), `tap_pu`, and `shift_deg`.
    Transformer tap values are in per unit on the listed bus voltage bases.
    """

    def __init__(self, base_mva: float, branch_parameters: pd.DataFrame):
        if base_mva <= 0:
            raise ValueError("base_mva must be positive")
        self.base_mva = float(base_mva)
        self.base_dir = ROOT
        self.bus_table = pd.read_csv(ROOT / "IEEE9BUS_loadflow_results.csv")
        self.topology = pd.read_csv(ROOT / "IEEE9_connections_loadflow.csv")
        required = {"name", *PARAMETER_COLUMNS}
        missing = required.difference(branch_parameters.columns)
        if missing:
            raise ValueError(f"branch parameter table lacks columns: {sorted(missing)}")
        params = branch_parameters.copy()
        if params["name"].duplicated().any():
            raise ValueError("branch parameter names must be unique")
        branches = self.topology.merge(params, on="name", how="left", validate="one_to_one")
        if branches[list(PARAMETER_COLUMNS)].isna().any().any():
            absent = branches.loc[branches[list(PARAMETER_COLUMNS)].isna().any(axis=1), "name"].tolist()
            raise ValueError(f"missing physical parameters for branches: {absent}")
        if set(params["name"]) != set(self.topology["name"]):
            raise ValueError("branch parameter names must match the 9 topology rows exactly")
        if (branches["r_pu"] < 0).any() or (branches["x_pu"] == 0).any():
            raise ValueError("branch resistance must be nonnegative and reactance nonzero")
        if (branches["tap_pu"] <= 0).any():
            raise ValueError("tap_pu values must be positive")
        self.branches = branches
        self.bus_index = {bus: i for i, bus in enumerate(BUS_IDS)}
        self.ybus = self._build_ybus()
        self.yload = self._load_admittance()

    def _build_ybus(self):
        ybus = np.zeros((len(BUS_IDS), len(BUS_IDS)), dtype=complex)
        for row in self.branches.itertuples(index=False):
            i = self.bus_index[int(row.from_bus)]
            j = self.bus_index[int(row.to_bus)]
            y_series = 1.0 / complex(float(row.r_pu), float(row.x_pu))
            y_shunt = 0.5j * float(row.b_pu)
            tap = float(row.tap_pu) * np.exp(1j * np.deg2rad(float(row.shift_deg)))
            ybus[i, i] += (y_series + y_shunt) / (abs(tap) ** 2)
            ybus[i, j] -= y_series / np.conj(tap)
            ybus[j, i] -= y_series / tap
            ybus[j, j] += y_series + y_shunt
        return ybus

    def _load_admittance(self):
        """Return the constant-impedance load shunts from the load-flow point."""
        loads = {
            # Bus assignments are recovered by matching each load's listed
            # voltage/angle to the corresponding Busbar row in the workbook.
            5: (127.54299132456107, 50.995615185477895),
            6: (50.930209430376706 + 91.674376974678196,
                0.0 + 30.558125658225752),
            8: (107.03588678199104, 37.442753157276869),
        }
        out = np.zeros(len(BUS_IDS), dtype=complex)
        buses = self.bus_table.loc[self.bus_table["BlockType"] == "Busbar"].copy()
        for bus, (p_mw, q_mvar) in loads.items():
            row = buses.loc[buses["Row"].str.startswith(f"Bus{bus} ")]
            if len(row) != 1:
                raise ValueError(f"could not identify load-flow voltage for Bus{bus}")
            vm = float(row.iloc[0]["ActualVoltageMagnitude_Pu"])
            s_pu = complex(p_mw, q_mvar) / self.base_mva
            out[self.bus_index[bus]] = np.conj(s_pu) / (vm * vm)
        return out

    def algebraic_current_residual(self, voltage_pu, generator_current_injections_pu):
        """Return g(z,x,u,p)=0 current balance for nine complex bus voltages.

        Generator current injections are positive into the network. They must
        be supplied by a coupled generator model; this method only implements
        the network and constant-impedance load side of the DAE.
        """
        voltage = np.asarray(voltage_pu, dtype=complex)
        current = np.asarray(generator_current_injections_pu, dtype=complex)
        if voltage.shape != (9,) or current.shape != (9,):
            raise ValueError("voltage and generator current vectors must have length 9")
        return (self.ybus @ voltage) + self.yload * voltage - current

    def loadflow_power_mismatch(self, voltage_pu, generation_mva):
        """Compute bus complex-power mismatch using constant-Z loads.

        `generation_mva` is a length-nine complex vector of P+jQ injections.
        This is a pointwise load-flow diagnostic, not a dynamic validation.
        """
        voltage = np.asarray(voltage_pu, dtype=complex)
        generation = np.asarray(generation_mva, dtype=complex)
        if voltage.shape != (9,) or generation.shape != (9,):
            raise ValueError("voltage and generation vectors must have length 9")
        s_network_pu = voltage * np.conj(self.ybus @ voltage)
        s_load_pu = self.yload.real * np.abs(voltage) ** 2 + 1j * (
            -self.yload.imag * np.abs(voltage) ** 2
        )
        s_gen_pu = generation / self.base_mva
        return (s_network_pu + s_load_pu - s_gen_pu) * self.base_mva


def load_branch_parameters(path: str | Path):
    """Read an explicit branch-parameter CSV supplied by the user/model source."""
    return pd.read_csv(path)
