import numpy as np
import matplotlib.pyplot as plt

times = [20, 40, 60]
safety_factors = ["03", "06", "09"]

rho_max = 0.1
v_max = 50.0

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for col, t in enumerate(times):
    csv_name = f"rhovx_exact_{t}.csv"
    data = np.loadtxt(csv_name, delimiter=",")

    ax = axes[col]

    ax.plot(data[:,0], data[:,1])
    ax.axvline(0, color="green", linestyle="--")
    ax.set_ylim(0, 0.11)

    ax.set_xlabel("Distance [ft]")
    if col == 0:
        ax.set_ylabel("Density [Cars/ft]")

    ax.set_title(f"Exact Density at t = {t} s")

fig.tight_layout()
plt.savefig("exact_collage.png")
plt.close()

fig, axes = plt.subplots(3, 3, figsize=(16, 14))

for row, S_f in enumerate(safety_factors):
    for col, t in enumerate(times):

        csv_name = f"rhovx_{t}_{S_f}.csv"
        data = np.loadtxt(csv_name, delimiter=",")
        data_exact = np.loadtxt(f"rhovx_exact_{t}.csv", delimiter=",")

        rho_exact_interp = np.interp(data[:,0], data_exact[:,0], data_exact[:,1])
        rmse = np.sqrt(np.mean((data[:,1] - rho_exact_interp)**2))

        if data[21, 1] == 0:
            rho_light = data[20, 1]
        else:
            rho_light = (data[21, 1] + data[20, 1]) / 2
        v_light = v_max * (1 - rho_light / rho_max)

        idx_first = int(np.max(np.where(data[:,1] > 0.01)))
        rho_first = data[idx_first, 1]
        v_first = v_max * (1 - rho_first / rho_max)

        idx_last = int(np.min(np.where(data[:,1] > 0.01)))
        rho_last = data[idx_last + 1, 1]
        v_last = v_max * (1 - rho_last / rho_max)

        Sf_value = int(S_f) / 10

        ax = axes[row, col]

        ax.plot(data[:,0], data[:,1], linewidth=2)
        ax.axvline(0, color="green", linestyle="--", linewidth=1.5)
        ax.set_ylim(0, 0.11)

        if row == 2:
            ax.set_xlabel("Distance [ft]", fontsize=14)

        if col == 0:
            ax.set_ylabel(f"S_f = {Sf_value}\nDensity [cars/ft]", fontsize=14)

        ax.set_title(f"t = {t} s", fontsize=16)

        ax.text(0.98, 0.95, f"v @ light = {v_light:.2f} ft/s",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=12, fontweight="bold")

        ax.text(0.98, 0.88, f"v @ 1st car = {v_first:.2f} ft/s",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=12, fontweight="bold")

        ax.text(0.98, 0.81, f"v @ 100th car = {v_last:.2f} ft/s",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=12, fontweight="bold")

        ax.text(0.98, 0.74, f"RMSE = {rmse:.5f} cars/ft",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=12, fontweight="bold")

fig.tight_layout()
plt.savefig("collage.png")
plt.close()
