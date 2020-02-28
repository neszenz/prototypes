import matplotlib.pyplot as plt
import numpy as np
import pathlib

from meshkD import SuperVertex, MeshkD, load_from_file
from util import *

def load_min_angles_list(input_dir):
    def min_angle_in_mesh(mesh):
        min_angle = float('inf')

        for face_mesh in mesh.face_meshes:
            vertices, _, triangles, _ = face_mesh
            for (i0, i1, i2) in triangles:
                p0 = vertices[i0].XYZ_vec3()
                p1 = vertices[i1].XYZ_vec3()
                p2 = vertices[i2].XYZ_vec3()

                curr_min_angle = calculate_min_angle_in_triangle(p0, p1, p2)

                if curr_min_angle < min_angle:
                    min_angle = curr_min_angle

        return np.rad2deg(min_angle)
    print('loading', input_dir)
    angles = []

    paths = pathlib.Path(input_dir).rglob('*'+MeshkD.FILE_EXTENSION)
    for path in sorted(paths):
        mesh = load_from_file(path)
        angles.append(min_angle_in_mesh(mesh))

    return angles

def smooth_angles_list(angles_list, window_size):
    smooth_angles_list = []

    window_size = max(3, window_size)
    if window_size % 2 == 0:
        window_size += 1

    radius = int(window_size / 2)
    n_angles = len(angles_list)

    for i in range(n_angles):
        if i == 0 or i == n_angles-1:
            smooth_angles_list.append(angles_list[i])
            continue

        angle_sum = 0

        j = i-radius
        while j < i+radius+1:
            if j < 0 or j > n_angles-1:
                angle_sum += angles_list[i]
            else:
                angle_sum += angles_list[j]
            j += 1

        angle_average = angle_sum / window_size
        smooth_angles_list.append(angle_average)

    return smooth_angles_list

SIZE = 10
SMOOTH_WINDOW_SIZE = 5

AREA_DIR = '../results/thesis_figures/convergence_graph/area'
CIRCUMRADIUS_DIR = '../results/thesis_figures/convergence_graph/circumradius'
DISTANCE_DIR = '../results/thesis_figures/convergence_graph/distance'

smooth_angles_list(np.arange(10), SMOOTH_WINDOW_SIZE)

angles_area = load_min_angles_list(AREA_DIR)
angles_circumradius = load_min_angles_list(CIRCUMRADIUS_DIR)
angles_distance = load_min_angles_list(DISTANCE_DIR)

# smoothing by averaging in SMOOTH_WINDOW_SIZE around each value
angles_area = smooth_angles_list(angles_area, SMOOTH_WINDOW_SIZE)
angles_circumradius = smooth_angles_list(angles_circumradius, SMOOTH_WINDOW_SIZE)
angles_distance = smooth_angles_list(angles_distance, SMOOTH_WINDOW_SIZE)

min_angle_area = min(angles_area)
min_angle_circumradius = min(angles_circumradius)
min_angle_distance = min(angles_distance)
MIN_ANGLE = min(min_angle_area, min_angle_circumradius, min_angle_distance)

n_iter_area = len(angles_area) + 1
n_iter_circumradius = len(angles_circumradius) + 1
n_iter_distance = len(angles_distance) + 1
MAX_ITER = max(n_iter_area, n_iter_circumradius, n_iter_distance) - 1

# collect endpoints for scatter plot
endpointsX = [n_iter_area-1, n_iter_circumradius-1, n_iter_distance-1]
endpointsY = [angles_area[-1], angles_circumradius[-1], angles_distance[-1]]

# insert min angle in all, because initial triangulation was not stored
angles_area.insert(0, MIN_ANGLE)
angles_circumradius.insert(0, MIN_ANGLE)
angles_distance.insert(0, MIN_ANGLE)

fig, ax = plt.subplots(1, 1, figsize=(SIZE, SIZE/1.33))

ax.set_prop_cycle(color=['#1f77b4', '#ff7f0e', '#2ca02c'])
ax.spines['top'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['left'].set_visible(True)

ax.set(xlabel='refinement iteration',
       ylabel='minimum inner angle',
       ylim=(0, 35),
       xlim=(0, MAX_ITER+5),
       title='minimum angle convergence (smoothed) for different priority metrics')

ax.set_xticks(range(0, MAX_ITER, 20))
ax.set_yticks(range(0, 35, 5))
ax.yaxis.set_major_formatter(plt.FuncFormatter('{:.0f}°'.format))

ax.plot(np.arange(n_iter_area), angles_area, label='area')
ax.plot(np.arange(n_iter_circumradius), angles_circumradius, label='circumradius')
ax.plot(np.arange(n_iter_distance), angles_distance, label='distance')
ax.scatter(endpointsX, endpointsY, c='black', label='termination')

ax.grid(True, ls='-.')
ax.legend(loc='upper left')
fig.tight_layout()

plt.show()
