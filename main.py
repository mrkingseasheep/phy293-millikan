import math

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

FPS = 10
PX_TO_M = 540 * 1000  # 540 px / mm * 1000 mm / m
A = 5 / 4
lambda_path = 67.5e-9
eta = 1.827e-5
g = 9.81
rho = 875.3
d = 6e-3


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


def find_lin_fit(points):
    time_s = points[0]
    position_m = points[1]
    position_err_m = np.ones(len(time_s)) * (1 / PX_TO_M)
    popt, pcov = curve_fit(
        line, time_s, position_m, sigma=position_err_m, absolute_sigma=True
    )
    m, b = popt
    m_err, b_err = np.sqrt(np.diag(pcov))
    chi_squared = 0

    for i in range(len(time_s)):
        chi_squared += ((position_m[i] - (m * time_s[i] + b)) / position_err_m[i]) ** 2

    red_chi_squared = chi_squared / (len(time_s) - 2)
    p_value = stats.chi2.sf(chi_squared, len(time_s) - 2)

    return abs(m), abs(m_err), red_chi_squared, p_value


def find_cunningham_adjusted_radius(v_d, v_d_err):
    r = (
        -A * lambda_path
        + math.sqrt((A * lambda_path) ** 2 + 18 * eta * abs(v_d) / g / rho)
    ) / 2

    r_err = abs(
        (-9 * eta * v_d_err)
        / (
            2
            * g
            * rho
            * math.sqrt((A * lambda_path) ** 2 + 18 * eta * abs(v_d) / g / rho)
        )
    )

    return r, r_err


def find_charge(filename, voltage_rise, voltage_rise_err):
    rise_points, fall_points = read_file(filename)
    rise_v, rise_v_err, rise_redchi2, rise_p = find_lin_fit(rise_points)
    fall_v, fall_v_err, fall_redchi2, fall_p = find_lin_fit(fall_points)

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
    find_charge("data/0.txt", 460, 5)


if __name__ == "__main__":
    main()
