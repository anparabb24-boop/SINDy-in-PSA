import numpy as np
import pandas as pd

# ==========================================
# 1. LOAD DATASETS (SINGLE CSV + TIME SERIES)
# ==========================================
base_dir = "/Users/Arya/Desktop/gcn_env/sindy_in_PSA"

# Load load-flow results CSV
loadflow_df = pd.read_csv(f"{base_dir}/IEEE9BUS_loadflow_results.csv")

# Load generator time-series baselines
g1_df = pd.read_csv(f"{base_dir}/healthy datasets/G1_H.csv")
g2_df = pd.read_csv(f"{base_dir}/healthy datasets/G2_H.csv")
g3_df = pd.read_csv(f"{base_dir}/healthy datasets/G3_H.csv")

print("--- Data Successfully Loaded ---")
print(f"Load-flow dataframe shape: {loadflow_df.shape}")

# ==========================================
# 2. EXTRACT STEADY-STATE INITIALIZATION
# ==========================================
# Extract generator initial conditions from G1, G2, G3 CSVs at t = 0
init_states = {}
for name, df in zip(["Gen1", "Gen2", "Gen3"], [g1_df, g2_df, g3_df]):
  init_states[name] = {
      "t0": df["t"].iloc[0],
      "fd_p": df["fd_p"].iloc[0],  # Field excitation voltage
      "Pm_pu": df["Pm_pu"].iloc[0],  # Mechanical power input (pu)
      "theta": df["theta"].iloc[0],  # Rotor angle (deg)
      "omega": df["omega"].iloc[0],  # Rotor speed (rad/s)
  }

print("\n[Generator 1 Initial Conditions from G1_H.csv]:")
for k, v in init_states["Gen1"].items():
  print(f"  {k}: {v}")

# Filter busbars and synchronous machines from the load-flow data
busbars_df = loadflow_df[loadflow_df["BlockType"] == "Busbar"]
generators_df = loadflow_df[loadflow_df["BlockType"] == "Synchronous Machine"]

print("\n[Bus Voltage Profiles (First 3 Buses)]:")
print(
    busbars_df[
        [
            "Row",
            "RatedVoltage_KV",
            "ActualVoltageMagnitude_Pu",
            "VoltageAngle_Deg",
        ]
    ].head(3)
)


# ==========================================
# 3. SETUP NETWORK TOPOLOGY & DAE TEMPLATE
# ==========================================
class IEEE9BusDAESystem:

  def __init__(self, loadflow_df, generator_inits):
    self.loadflow = loadflow_df
    self.generator_inits = generator_inits
    self.system_frequency = 60.0  # Hz
    self.omega_s = 2 * np.pi * self.system_frequency  # Synchronous speed (377 rad/s)

  def initialize_state_vector(self):
    """Constructs the initial DAE state vector x0 combining generators and network nodes."""
    x0 = []
    for gen_name in ["Gen1", "Gen2", "Gen3"]:
      delta_0 = np.radians(self.generator_inits[gen_name]["theta"])
      omega_0 = self.generator_inits[gen_name]["omega"]
      fd_0 = self.generator_inits[gen_name]["fd_p"]

      # Append generator states [rotor angle, rotor speed, field voltage]
      x0.extend([delta_0, omega_0, fd_0])
    return np.array(x0)

  def compute_algebraic_residuals(self, t, x, V_bus, theta_bus):
    """Placeholder for network algebraic power balance equations (Y_bus coupling)."""
    residuals = []
    return residuals


# Instantiate and prepare system
dae_system = IEEE9BusDAESystem(loadflow_df, init_states)
x0 = dae_system.initialize_state_vector()

print(f"\nInitial State Vector x0 shape: {x0.shape}")
print("Initialized DAE Base Successfully.")