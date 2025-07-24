import numpy as np
from pathlib import Path
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

output_folder = Path("lensing_table")
output_folder.mkdir(exist_ok=True, parents=True)

u_min = 0.3
u_max = 0.5  # < 0.568 to ignore multiple images
N_u_grid = 128

cos_psi_min = -0.9
cos_psi_max = 1.0
N_cos_psi_grid = 512

d_alpha = np.pi / (N_cos_psi_grid * 8)

cos_alpha_min = -0.2
cos_alpha_max = 1.0
N_cos_alpha_grid = 512

output_u_grid = np.linspace(u_min, u_max, N_u_grid)
with open(output_folder / "u.txt", "w") as out_file:
    out_file.write(f"{u_min} {u_max} {N_u_grid}\n")

output_cos_psi_grid = np.linspace(cos_psi_min, cos_psi_max, N_cos_psi_grid)
with open(output_folder / "cos_psi.txt", "w") as out_file:
    out_file.write(f"{cos_psi_min} {cos_psi_max} {N_cos_psi_grid}\n")

output_cos_alpha_grid = np.linspace(cos_alpha_min, cos_alpha_max, N_cos_alpha_grid)
with open(output_folder / "cos_alpha.txt", "w") as out_file:
    out_file.write(f"{cos_alpha_min} {cos_alpha_max} {N_cos_alpha_grid}\n")


def f(x, u, alpha):
    s2 = np.sin(alpha) ** 2
    c2 = 1 - s2
    x2 = x**2
    q = (2 - x2 - u * (1 - x2) ** 2 / (1 - u)) * s2
    tmp1 = np.sqrt(c2 + x2 * q)
    ret1 = x / tmp1
    return ret1


def cal_psi1(u, alpha):
    res, err = quad(f, 0, 1, args=(u, alpha))
    s = np.sin(alpha)
    ret1 = 2 * s / np.sqrt(1 - u) * res
    return ret1


def cal_psi2(u, alpha):
    s = np.sin(alpha)
    tmp = (2 * s) / np.sqrt(3 * (1 - u))
    p_over_R = -tmp * np.cos((np.arccos(3 * u / tmp) + 2 * np.pi) / 3)
    psi_max = cal_psi1(u / p_over_R, np.pi / 2)
    psi1 = cal_psi1(u, np.pi - alpha)
    return 2 * psi_max - psi1


def cal_psi(u, alpha):
    if alpha <= np.pi / 2:
        return cal_psi1(u, alpha)
    else:
        return cal_psi2(u, alpha)


output_psi_grid = np.arccos(output_cos_psi_grid)
cos_alpha_of_u_cos_psi = np.zeros((N_u_grid, N_cos_psi_grid))
lf_of_u_cos_psi = np.zeros((N_u_grid, N_cos_psi_grid))

for u_id, u in enumerate(output_u_grid):
    print(f"Processing u = {u:.3f} ({u_id + 1}/{len(output_u_grid)})")
    alpha_grid, psi_grid = [], []
    alpha = 0
    while True:
        alpha_grid.append(alpha)
        psi_grid.append(cal_psi(u, alpha))
        alpha += d_alpha
        if psi_grid[-1] > output_psi_grid[0]:
            break
    alpha_grid = np.array(alpha_grid)
    psi_grid = np.array(psi_grid)
    cos_alpha_grid = np.cos(alpha_grid)
    cos_psi_grid = np.cos(psi_grid)
    # print(f"{cos_alpha_grid[-1]=}")

    lensing_factor = np.gradient(cos_alpha_grid, cos_psi_grid)

    cos_alpha_of_cos_psi = PchipInterpolator(cos_psi_grid[::-1], cos_alpha_grid[::-1])
    lf_of_cos_psi = PchipInterpolator(cos_psi_grid[::-1], lensing_factor[::-1])

    cos_alpha_of_u_cos_psi[u_id, :] = cos_alpha_of_cos_psi(output_cos_psi_grid)
    lf_of_u_cos_psi[u_id, :] = lf_of_cos_psi(output_cos_psi_grid)


np.savetxt(output_folder / "cos_alpha_of_u_cos_psi.txt", cos_alpha_of_u_cos_psi)
np.savetxt(output_folder / "lf_of_u_cos_psi.txt", lf_of_u_cos_psi)

#######################################################################################


def g(x, u, alpha):
    s2 = np.sin(alpha) ** 2
    c2 = 1 - s2
    x2 = x**2
    q = (2 - x2 - u * (1 - x2) ** 2 / (1 - u)) * s2
    tmp1 = np.sqrt(c2 + x2 * q)
    ret2 = x / (tmp1 * (1 + tmp1))
    return ret2


def cal_cdt_over_R1(u, alpha):
    res, err = quad(g, 0, 1, args=(u, alpha))
    s = np.sin(alpha)
    ret2 = 2 * s**2 / (1 - u) * res
    return ret2


def cal_cdt_over_R2(u, alpha):
    s = np.sin(alpha)
    tmp = (2 * s) / np.sqrt(3 * (1 - u))
    p_over_R = -tmp * np.cos((np.arccos(3 * u / tmp) + 2 * np.pi) / 3)
    dt_max = cal_cdt_over_R1(u / p_over_R, np.pi / 2)
    dt1 = cal_cdt_over_R1(u, np.pi - alpha)
    return 2 * dt_max - dt1


def cal_cdt_over_R(u, alpha):
    if alpha <= np.pi / 2:
        return cal_cdt_over_R1(u, alpha)
    else:
        return cal_cdt_over_R2(u, alpha)


cdt_over_R_of_u_cos_alpha = np.zeros((N_u_grid, N_cos_alpha_grid))
for u_id, u in enumerate(output_u_grid):
    print(f"Processing u = {u:.3f} ({u_id + 1}/{len(output_u_grid)})")
    for alpha_id, cos_alpha in enumerate(output_cos_alpha_grid):
        cdt_over_R_of_u_cos_alpha[u_id, alpha_id] = cal_cdt_over_R(
            u, np.arccos(cos_alpha)
        )
np.savetxt(output_folder / "cdt_over_R_of_u_cos_alpha.txt", cdt_over_R_of_u_cos_alpha)
