#!/usr/bin/python
from random import random

FRAME_NUM = 10
OBJ_NUM = 1024

def generate_matrix(offset):
    matrix = []
    for i in range(4):
        for j in range(4):
            if i == j:
                matrix.append(1.0)
            elif j == 3:
                matrix.append(offset[i])
            else:
                matrix.append(0.0)
    return matrix

def move_alone_x(stride):
    return generate_matrix([stride, 0, 0])

def move_alone_y(stride):
    return generate_matrix([0, stride, 0])

def move_alone_z(stride):
    return generate_matrix([0, 0, stride])

def rand_move():
    return generate_matrix([random(), random(), random()])

def print_matrix(matrix, fp):
    for dp in matrix:
        fp.write(str(dp) + "\n")


def two_batch():
    description = "\n\n\
This case separate all the object into two parts, both move over y axis.\n\
The two batch of objects collide once with each other in the middle.\n\n"
    print(description)

    op = open("y2batch.inp", "w")

    for i in range(FRAME_NUM):
        # half move downward
        for j in range(OBJ_NUM/2):
            print_matrix(move_alone_y(float(-0.25 * i)), op)
        # half move upward
        for j in range(OBJ_NUM/2):
            print_matrix(move_alone_y(float(0.25 * i)), op)

def rand_mov():
    description = "\n\n\
This case generates random movements.\n\n"
    print(description)

    op = open("rand.inp", "w")

    for i in range(FRAME_NUM):
        # half move downward
        for j in range(OBJ_NUM/2):
            print_matrix(move_alone_y(float(-0.25 * i)), op)
        # half move upward
        for j in range(OBJ_NUM/2):
            print_matrix(move_alone_y(float(0.25 * i)), op)

# two_batch()
rand_mov()