import math

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

OUTPUT_FILE_PATH = "results.txt"
FPS = 10
PX_TO_M = 540 * 1000  # 540 px / mm * 1000 mm / m
VOLTAGE_ERR = 5

START_IDX = 0
END_IDX = 1

A = 5 / 4
lambda_path = 67.5e-9
eta = 1.827e-5
g = 9.81
rho = 875.3
d = 6e-3


def read_voltage():
    rise_voltage, fall_voltage = [], []
    with open("data/voltage.txt", "r") as file:
        for i, line in enumerate(file):
            voltages = line.strip().split("\t")
            fall_voltage.append(int(float(voltages[0])))
            rise_voltage.append(int(float(voltages[1])))
    return np.array(fall_voltage), np.array(rise_voltage)


def read_file(filename):
    rise_t_s = ()
    fall_t_s = ()
    rise_pos_m = []
    rise_time_s = []
    fall_pos_m = []
    fall_time_s = []
    with open(filename, "r") as file:
        for i, line in enumerate(file):
            line = line.strip()
            if i == 0:
                line = line.split(" ")
                rise_t_s = (int(float(line[0])) / FPS, int(float(line[1])) / FPS)
                fall_t_s = (int(float(line[2])) / FPS, int(float(line[3])) / FPS)
                continue
            frame, px = line.split(" ")
            time_s = int(frame) / FPS
            distance_m = int(px) / PX_TO_M
            # do we do anything related to this? idk ;-;
            if time_s > rise_t_s[0] and time_s < rise_t_s[1]:
                rise_time_s.append(time_s)
                rise_pos_m.append(distance_m)
            if time_s > fall_t_s[0] and time_s < fall_t_s[1]:
                fall_time_s.append(time_s)
                fall_pos_m.append(distance_m)
    return (np.array(rise_time_s), np.array(rise_pos_m)), (
        np.array(fall_time_s),
        np.array(fall_pos_m),
    )


def line(x, m, b):
    return m * x + b


def plot_residuals(time_s, position_m, m, b):
    residual = position_m - (m * time_s + b)
    # for i in range(len(time_s)):
    #     residual.append(position_m - (m * time_s + b))
    # residual = np.array(residual)
    # time_s = np.array(time_s)
    fig, ax = plt.subplots()
    ax.bar(time_s, residual, width=0.1)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Residual (m)")
    # plt.savefig("residual_fall_0.png", dpi=300, bbox_inches="tight")
    plt.show()


def find_lin_fit(points):
    time_s = points[0]
    position_m = points[1]
    # calculate pos error with 1 pixel as absolute error
    position_err_m = np.ones(len(time_s)) * (1 / PX_TO_M)
    # calculate pos error as percentage uncertainty converted to abs error
    # for i, position in enumerate(position_m):
    #     position_err_m[i] = position * 1 / 540
    popt, pcov = curve_fit(
        line, time_s, position_m, sigma=position_err_m, absolute_sigma=True
    )
    m, b = popt
    m_err, b_err = np.sqrt(np.diag(pcov))
    chi_squared = 0

    residuals = position_m - line(time_s, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((position_m - np.mean(position_m)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    rmse = np.sqrt(np.mean(residuals**2))
    print("r^2")
    print(r_squared)
    print("rmse")
    print(rmse)
    print("m")
    print(m)
    print("actual err")
    print(1 / PX_TO_M)

    for i in range(len(time_s)):
        chi_squared += ((position_m[i] - (m * time_s[i] + b)) / position_err_m[i]) ** 2

    red_chi_squared = chi_squared / (len(time_s) - 2)
    p_value = stats.chi2.sf(chi_squared, len(time_s) - 2)

    # plot_residuals(time_s, position_m, m, b)

    # brownian motion adjustment
    delta_t = abs(time_s[0] - time_s[-1])
    k_b = 1.380649e-23
    T = 20 + 273.15
    r = 6.60141439066232e-07  # back calculated: average of all drop particle radius from no brownian motion
    brownian_correction_dist = math.sqrt(2 * k_b * T / 6 / np.pi / eta / r * delta_t)
    # CHANGE ME TO ADD BROWNIAN CORRECTION
    m = (m * delta_t + brownian_correction_dist) / delta_t

    return abs(m), abs(m_err), red_chi_squared, p_value


def find_cunningham_adjusted_radius(v_d, v_d_err):
    r = (
        -A * lambda_path
        + math.sqrt((A * lambda_path) ** 2 + 18 * eta * abs(v_d) / g / rho)
    ) / 2
    # print(r)

    r_err = abs(
        (-9 * eta * v_d_err)
        / (
            2
            * g
            * rho
            * math.sqrt((A * lambda_path) ** 2 + 18 * eta * abs(v_d) / g / rho)
        )
    )
    # print(r_err)
    return r, r_err


def find_charge(filename, voltage_rise, voltage_rise_err):
    rise_points, fall_points = read_file(filename)
    rise_v, rise_v_err, rise_redchi2, rise_p = find_lin_fit(rise_points)
    fall_v, fall_v_err, fall_redchi2, fall_p = find_lin_fit(fall_points)
    # print(f"{rise_redchi2}")
    # print(f"{filename} rise velocity has reduced chi2 of {rise_redchi2}")

    r, r_err = find_cunningham_adjusted_radius(fall_v, fall_v_err)

    C = 4 / 3 * np.pi * rho * g * d
    q = C * r**3 * (fall_v + rise_v) / voltage_rise / fall_v

    partial_q_r = C * 3 * r**2 * (fall_v + rise_v) / voltage_rise / fall_v
    partial_q_V_rise = -C * r**3 * (fall_v + rise_v) / voltage_rise**2 / fall_v
    partial_q_v_u = C * r**3 / voltage_rise / fall_v
    partial_q_v_d = -C * rise_v * r**3 / voltage_rise / fall_v**2
    # print(f"q / r       {partial_q_r}")
    # print(f"q / V       {partial_q_V_rise}")
    # print(f"q / v_u     {partial_q_v_u}")
    # print(f"q / v_d     {partial_q_v_d}")

    q_err = math.sqrt(
        (partial_q_r * r_err) ** 2
        + (partial_q_v_d * fall_v_err) ** 2
        + (partial_q_v_u * rise_v_err) ** 2
        + (partial_q_V_rise * voltage_rise_err) ** 2
    )

    # print(q)
    # print(q_err)
    return q, q_err


def main():
    fall_voltage, rise_voltage = read_voltage()
    with open(OUTPUT_FILE_PATH, "w") as file:
        file.write(f"q q_err\n")
        for i in range(START_IDX, END_IDX):
            filepath = f"data/{i}.txt"
            q, q_err = find_charge(filepath, rise_voltage[i], VOLTAGE_ERR)
            file.write(f"{q} {q_err}\n")


if __name__ == "__main__":
    main()
