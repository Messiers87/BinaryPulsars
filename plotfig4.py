import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("data/fig4_data.csv")
Po = sorted(df["Po_day"].unique())
Pp = sorted(df["Pp_s"].unique())
Z = df.pivot(index="Pp_s", columns="Po_day", values="gamma1").values

plt.figure(figsize=(6,5))
CS = plt.contourf(Po, Pp, Z, levels=20, cmap="viridis")
plt.colorbar(CS, label="gamma1_m")
plt.xlabel("Orbital period Po [days]")
plt.ylabel("Spin period Pp [s]")
plt.yscale("log")
plt.tight_layout()
plt.savefig("plots/fig4a_thirdTry.png", dpi=300)
plt.show()


