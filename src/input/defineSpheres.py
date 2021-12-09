#!/usr/bin/python
from random import random

OBJ_NUM = 1024
RADIUS = 0.045
STRIDE = 0.1


def y2batch(op):
    y = 0.5
    z = 0
    for i in range(OBJ_NUM/2):
        # x range = [-(OBJ_NUM/4) * stride, (OBJ_NUM/4) * stride]
        benchmark = OBJ_NUM/4
        x = (i - benchmark) * STRIDE
        op.write(str(x) + " " + str(y) + " " + str(z) + " " + str(RADIUS) + "\n")
    for i in range(OBJ_NUM/2):
        # x range = [-(OBJ_NUM/4) * stride, (OBJ_NUM/4) * stride]
        benchmark = OBJ_NUM/4
        x = (i - benchmark) * STRIDE
        op.write(str(x) + " " + str(-y) + " " + str(z) + " " + str(RADIUS) + "\n")


def randpos(op):
    for i in range(OBJ_NUM):
        for j in range(3):
            sign = random() > 0.5
            val = random() * 2
            dp = val
            if sign:
                dp = -dp
            op.write(str(dp) + " ")
        rad = random()
        op.write(str(rad) + "\n")

op = open("spheres.txt", "w")

y2batch(op)
# randpos(op)
