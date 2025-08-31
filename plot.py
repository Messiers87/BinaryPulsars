import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# File names
file_fig2 = "data/fig2_velocity.csv"
file_fig3 = "data/fig3_acc.csv"

# ---------------------------
# Safety checks
# ---------------------------
missing_files = []
for fname in [file_fig2, file_fig3]:
    if not os.path.exists(fname):
        missing_files.append(fname)

if missing_files:
    print(" Error: The following required CSV files were not found:")
    for f in missing_files:
        print("   -", f)
    print("\nGenerate plots if not already by running c++ program")
    sys.exit(1)

print(" Found required CSV files, proceeding with plotting...")

# ---------------------------
# Figure 2: velocity vs time
# ---------------------------
df2 = pd.read_csv(file_fig2)

plt.figure(figsize=(10,8))

# Upper panel
plt.subplot(2,1,1)
plt.plot(df2["t_sec"], df2["v_varpi_0deg"], label="varpi=0°")
plt.plot(df2["t_sec"], df2["v_varpi_90deg"], label="varpi=90°")
plt.plot(df2["t_sec"], df2["v_varpi_180deg"], label="varpi=180°")
plt.plot(df2["t_sec"], df2["v_varpi_270deg"], label="varpi=270°")
plt.plot(df2["t_sec"], df2["v_circular"], label="circular", linestyle="--")
plt.xlabel("t (s)")
plt.ylabel("v_l (m/s)")
plt.title("Figure 2 (Upper panel)")
plt.legend()

# Lower panel
plt.subplot(2,1,2)
plt.plot(df2["t_sec"], df2["v_varpi_60deg"], label="varpi=60°")
plt.plot(df2["t_sec"], df2["v_varpi_120deg"], label="varpi=120°")
plt.plot(df2["t_sec"], df2["v_varpi_240deg"], label="varpi=240°")
plt.plot(df2["t_sec"], df2["v_varpi_300deg"], label="varpi=300°")
plt.plot(df2["t_sec"], df2["v_circular"], label="circular", linestyle="--")
plt.xlabel("t (s)")
plt.ylabel("v_l (m/s)")
plt.title("Figure 2 (Lower panel)")
plt.legend()

plt.tight_layout()
plt.savefig("plots/fig2_velocity_py.png", dpi=200)
plt.close()

# ---------------------------
# Figure 3: acceleration vs f
# ---------------------------
df3 = pd.read_csv(file_fig3)

plt.figure(figsize=(8,5))
plt.plot(df3["f_deg"], df3["a_l__m/s2"], label="a_l(f)")
plt.xlabel("True anomaly f (deg)")
plt.ylabel("a_l (m/s²)")
plt.title("Figure 3")
plt.legend()
plt.grid(True)

plt.savefig("plots/fig3_acc_py.png", dpi=200)
plt.close()

print(" Saved fig2_vel_py.png and fig3_accel_py.png")
