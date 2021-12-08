#!/usr/bin/python
import math

PI = 3.14159
THRESHOLD = 0.0001
SLICE = 4
STACK = 8

def my_round(f):
    return f
    # rounded_f = round(f)
    # if abs(rounded_f - f) < THRESHOLD:
    #     return rounded_f
    # else:
    #     return f

def degree_to_radian(d):
    return PI * d / 180

# record the upper part of the sphere as points
def sphere_to_point(org, rad, gratitude, lattitude):
    round_buffer = []

    for d_phi in range(0, 90, 90/gratitude):
        phi = degree_to_radian(90 - d_phi)
        upper_buffer = []
        for d_theta in range(0, 360, 360/lattitude):
            theta = degree_to_radian(d_theta)
            p = []
            x = my_round(org[0] + rad * math.sin(phi) * math.cos(theta))
            p.append(x)
            y = my_round(org[1] + rad * math.sin(phi) * math.sin(theta))
            p.append(y)
            z = my_round(org[2] + rad * math.cos(phi))
            p.append(z)

            # print(x, y, z)

            upper_buffer.append(p)
        upper_buffer.append(upper_buffer[0])    # close the circle
        round_buffer.insert(0, upper_buffer)

    return round_buffer

def mirrowP(point, org_z):
    # print("old:", point)
    mirrow_point = []
    for  i in range(3):
        mirrow_point.append(point[i])
    if point[2] != 0:
        mirrow_point[2] = org_z-point[2]
    # print("new:", mirrow_point)

    return mirrow_point

def mirrowT(triangle, org_z):
    new_triangle = []
    # print("new triangle ----")
    for p in triangle:
        new_p = mirrowP(p, org_z)
        new_triangle.append(new_p)
    return new_triangle

def partialMirrowT(triangle, org_z):
    new_triangle = []
    # print("new triangle ----")
    for i in range(2):
        p = triangle[i]
        new_p = mirrowP(p, org_z)
        new_triangle.append(new_p)
    new_triangle.append(p)
    return new_triangle

def point_to_triangle(round_buffer, top, org_z):
    upper_buffer = []
    lower_buffer = []

    # deal with the highest gratitude -> render triangles
    cur_round_upper = []  # list of triangles of the current_round for the upper half of the sphere
    cur_round_lower = []  # list of triangles of the current_round for the bottom half of the sphere
    cur_round = round_buffer[0]
    for i in range(1, len(cur_round)):
        p_left = cur_round[i - 1]
        p_right = cur_round[i]
        triangle = [top, p_left, p_right]

        # print("tri:", triangle)
        # print("mir:", mirrowT(triangle, org_z))

        cur_round_upper.append(triangle)
        mirrow = mirrowT(triangle, org_z)
        cur_round_lower.append(mirrow)
        # print("triangle", triangle, mirrow)
        # print("upper before", cur_round_upper)
        # print("lower before", cur_round_lower)
        
    upper_buffer.append(cur_round_upper)
    lower_buffer.insert(0, cur_round_lower)

    # deal with the lower gratitudes -> render rectangles -> separate to triangles
    for j in range(1, len(round_buffer)):  # j -> gratitude
        edge = (len(round_buffer) - 1 == j)

        round_high = round_buffer[j - 1]  # point
        round_low = round_buffer[j]
        cur_round_upper = []            # triangle
        cur_round_lower = []
        for i in range(len(cur_round)):
            i_left = i
            i_right = (i + 1) % len(cur_round)

            top_left = round_high[i_left]
            top_right = round_high[i_right]
            bottom_left = round_low[i_left]
            bottom_right = round_high[i_right]

            tri_left = [top_left, top_right, bottom_left]
            tri_right = [top_right, bottom_left, bottom_right]
            cur_round_upper.append(tri_left)
            cur_round_upper.append(tri_right)

            if edge:
                cur_round_lower.append(partialMirrowT(tri_left, org_z))
                cur_round_lower.append(partialMirrowT(tri_right, org_z))
            else:
                cur_round_lower.append(mirrowT(tri_left, org_z))
                cur_round_lower.append(mirrowT(tri_right, org_z))

        upper_buffer.append(cur_round_upper)
        lower_buffer.insert(0, cur_round_lower)

    triangle_buffer = upper_buffer + lower_buffer

    # flatten the records
    triangle_list = []
    for round_buffer in triangle_buffer:
        for tri in round_buffer:
            triangle_list.append(tri)
    
    return triangle_list


# extract sphere info
lines = []
with open('spheres.txt') as f:
    lines = f.readlines()

tri_count = 0
obj_buffer = []
for line in lines:
    dim = line.split(" ")
    org = [float(dim[0]), float(dim[1]), float(dim[2])]
    rad = float(dim[3])
    
    # process to find all triangles for each sphere
    point_list = sphere_to_point(org, rad, SLICE, STACK)
    triangle_list = point_to_triangle(point_list, [org[0], org[1], org[2] + rad], org[2])
    tri_count = len(triangle_list)     # triangle count for one sphere
    obj_buffer.append(triangle_list)

# write to create inp file
op = open("triangle.inp", "w")
op.write(str(tri_count))
op.write('\n\n')

for obj in obj_buffer:
    for tri in obj:
        for p in tri:
            op.write(str(p[0]) + " " + str(p[1]) + " " + str(p[2]) + "\n")
        op.write("\n")
 