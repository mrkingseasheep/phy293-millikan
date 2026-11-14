import math

import matplotlib.pyplot as plt
import numpy as np

FPS = 10
PX_TO_M = 540 * 1000  # 540 px / mm * 1000 mm / m


def read_file(filename):
    rise_t_s = ()
    fall_t_s = ()
    rise_pos_s_m = []
    fall_pos_s_m = []
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
                rise_pos_s_m.append((time_s, distance_m))
            if time_s > fall_t_s[0] and time_s < fall_t_s[1]:
                fall_pos_s_m.append((time_s, distance_m))
    return rise_pos_s_m, fall_pos_s_m


def main():
    rise_pos_s_m, fall_pos_s_m = read_file("data/0.txt")
    print(rise_pos_s_m)
    print(fall_pos_s_m)


main()
