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

def calculate_cog_2d(sv0, sv1, sv2):
    return (sv0.UV_vec2() + sv1.UV_vec2() + sv2.UV_vec2()) / 3

def calculate_cog_3d(sv0, sv1, sv2):
    return (sv0.XYZ_vec3() + sv1.XYZ_vec3() + sv2.XYZ_vec3()) / 3

def calculate_cog(sv0, sv1, sv2):
    return calculate_cog_3d(sv0, sv1, sv2), calculate_cog_2d(sv0, sv1, sv2)

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

def plot_histogram(values, name, num_of_bins, bounds):
    lower_bound, upper_bound = bounds
    if lower_bound == float('-inf'):
        bounds = (min(values), upper_bound)
    if upper_bound == float('inf'):
        bounds = (lower_bound, max(values))

    n, bins, patches = plt.hist(values, num_of_bins, bounds, facecolor='blue', alpha=0.5)
    plt.suptitle(name)
    plt.show()

    return

# custom metric for measuring triangle regularity via normalized area
def calculate_xi(p0, p1, p2):
    l0 = np.linalg.norm(p1 - p0)
    l1 = np.linalg.norm(p2 - p1)
    l2 = np.linalg.norm(p0 - p2)
    assert min(l0, l1, l2) > 0.0

    # normalize for longest edge
    lmax = max(l0, l1, l2)
    l0 /= lmax
    l1 /= lmax
    l2 /= lmax

    # xi now is the area (via heron's formula) for normalized side lengths
    s = 0.5 * (l0 + l1 + l2)
    xi = np.sqrt(s*(s-l0)*(s-l1)*(s-l2))

    return xi

def calculate_circumradius_v1(p0, p1, p2):
    # based on book 'Delaunay Mesh Generation' page 26
    l01 = np.linalg.norm(p1 - p0)
    l12 = np.linalg.norm(p2 - p1)
    l20 = np.linalg.norm(p0 - p2)

    l_min = l01 if l01 < l12 else l12
    l_min = l20 if l20 < l_min else l_min

    alpha = calculate_angle_in_corner(p2, p0, p1)
    beta = calculate_angle_in_corner(p0, p1, p2)
    gamma = calculate_angle_in_corner(p1, p2, p0)

    a_min = alpha if alpha < beta else beta
    a_min = gamma if gamma < a_min else a_min
    sin_a_min = np.sin(a_min)
    assert sin_a_min > 0.0

    return l_min / (2 * np.sin(a_min))

def calculate_circumradius_v2(p0, p1, p2):
    # based on https://artofproblemsolving.com/wiki/index.php/Circumradius
    l01 = np.linalg.norm(p1 - p0)
    l12 = np.linalg.norm(p2 - p1)
    l20 = np.linalg.norm(p0 - p2)

    p01 = p1 - p0
    p02 = p2 - p0
    double_area = np.linalg.norm(np.cross(p01, p02))
    assert double_area > 0.0

    return (l01*l12*l20) / (2*double_area)

def calculate_radius_edge_ratio_v1(p0, p1, p2):
    # based on https://artofproblemsolving.com/wiki/index.php/Circumradius
    l01 = np.linalg.norm(p1 - p0)
    l12 = np.linalg.norm(p2 - p1)
    l20 = np.linalg.norm(p0 - p2)

    p01 = p1 - p0
    p02 = p2 - p0
    double_area = np.linalg.norm(np.cross(p01, p02))
    assert double_area > 0.0

    circumradius = (l01*l12*l20) / (2*double_area)
    l_min = min(l01, l12, l20)
    assert l_min > 0.0

    return circumradius / l_min

def calculate_radius_edge_ratio_v2(p0, p1, p2):
    # based on 'Delaunay Mesh Generation' page 26
    a = calculate_angle_in_corner(p2, p0, p1)
    b = calculate_angle_in_corner(p0, p1, p2)
    c = calculate_angle_in_corner(p1, p2, p0)

    theta_min = min(a, b, c)

    return 1.0 / (2.0 * np.sin(theta_min))
