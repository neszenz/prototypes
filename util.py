import numpy as np
import matplotlib.pyplot as plt

from OCC.Core.BRepAdaptor import BRepAdaptor_Surface

def normalize(vector):
    vNorm = np.linalg.norm(vector)
    if vNorm == 0.0:
        return vector
    else:
        return vector / vNorm

def calculate_angle_between_vectors(v0, v1):
    dp = max(-1.0, min(1.0, np.dot(v0, v1)))

    return np.arccos(dp)

def calculate_angle_in_corner(p0, p1, p2):
    p10 = normalize(p0 - p1)
    p12 = normalize(p2 - p1)

    return calculate_angle_between_vectors(p10, p12)

def calculate_normal(p0, p1, p2):
    p10 = p0 - p1
    p12 = p2 - p1

    return np.cross(p12, p10)

def calculate_normal_normalized(p0, p1, p2):
    return normalize(calculate_normal(p0, p1, p2))

def calculate_area(p0, p1, p2):
    return np.linalg.norm(calculate_normal(p0, p1, p2)) / 2.0

def left_hand_perpendicular(v):
    assert len(v) == 2
    return np.array((v[1], -v[0]))

def project_point_onto_normalized_vector(p, v):
    return np.dot(p, v) * v

# assuming all arguments are numpy arrays and the directions are normalized
def shortest_vector_between_two_lines(line0_ori, line0_dir, line1_ori, line1_dir):
    # set up and solve linear equation
    closest_connection_dir = np.cross(line0_dir, line1_dir)
    A = np.array([line0_dir, -line1_dir, closest_connection_dir]).T
    B = line1_ori - line0_ori
    x = np.linalg.solve(A, B)

    # calculate closest points
    line0_fact, line1_fact, connection_fact = x
    p0 = line0_ori + (line0_fact * line0_dir)
    p1 = line1_ori + (line1_fact * line1_dir)

    return p1 - p0

def reverse_u(u, face):
    last_u_parameter = BRepAdaptor_Surface(face).LastUParameter()
    return last_u_parameter - u

def plot_histogram(values, num_of_bins, lower_bound=None, upper_bound=None):
    if lower_bound is None or upper_bound is None:
        n, bins, patches = plt.hist(values, num_of_bins, facecolor='blue', alpha=0.5)
    else:
        n, bins, patches = plt.hist(values, num_of_bins, range=(lower_bound, upper_bound), facecolor='blue', alpha=0.5)
    plt.show()
