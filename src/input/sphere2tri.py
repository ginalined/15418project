#!/usr/bin/python
import math

M_PI = 3.14159

def transfer(org, rad, n_stack, n_slice):
    p_buffer = []
    for i in range(n_stack):
        phi = M_PI * (i + 1) / (n_stack)
        for j in range(n_slice):
            p = []
            theta = 2 * rad * M_PI * j / n_slice

            x = org[0] + rad * math.sin(phi) * math.cos(theta)
            p.append(x)
            y = org[1] + rad * math.sin(phi) * math.sin(theta)
            p.append(y)
            z = org[2] + rad * math.cos(phi)
            p.append(z)

            print("p: ", p)
            p_buffer.append(p)
    
    triangle_buffer = []
    
    v0 = [org[0], org[1], org[2] + rad] # top
    p_buffer.append(v0)
    print("top: ", v0)
    v1 = [org[0], org[1], org[2] - rad] # bottom
    p_buffer.append(v1)
    print("bottom: ", v1)

    print(len(p_buffer))

    # add top triangles
    for i in range(n_slice):
        triangle = [v0]
        i0 = i + 1
        triangle.append(p_buffer[i0])
        i1 = (i + 1) % n_slice + 1
        triangle.append(p_buffer[i1])
        triangle_buffer.append(triangle)    

    # add quads per stack / slice
    for j in range(n_stack - 2):
        j0 = j * n_slice + 1
        j1 = (j + 1) * n_slice + 1
        for i in range(n_slice):
            i0 = j0 + i
            i1 = j0 + (i + 1) % n_slice
            i2 = j1 + (i + 1) % n_slice
            i3 = j1 + i

            print(i0, i1, i2, i3)

            tri_left = [p_buffer[i0], p_buffer[i1], p_buffer[i3]]
            tri_right = [p_buffer[i1], p_buffer[i3], p_buffer[i2]]

            triangle_buffer.append(tri_left)
            triangle_buffer.append(tri_right)

    # add bottom triangles
    for i in range(n_slice):
        triangle = [v1]
        i0 = i + n_slice * (n_stack - 2) + 1
        triangle.append(p_buffer[i0])
        i1 = (i + 1) % n_slice + n_slice * (n_stack - 2) + 1
        triangle.append(p_buffer[i1])
        triangle_buffer.append(triangle)
    
    return triangle_buffer



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
    triangle_list = transfer(org, rad, 3, 3)
    tri_count += len(triangle_list)
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
 