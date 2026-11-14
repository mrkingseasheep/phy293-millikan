import math

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

FPS = 10
PX_TO_M = 540 * 1000  # 540 px / mm * 1000 mm / m


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


def find_lin_fit(points):
    time_s = points[0]
    position_m = points[1]
    print(position_m)
    print(time_s)


def main():
    rise_points, fall_points = read_file("data/0.txt")
    # print(rise_pos_s_m)
    find_lin_fit(rise_points)


if __name__ == "__main__":
    main()
