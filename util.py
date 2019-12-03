import numpy as np

def normalize(vector):
    vNorm = np.linalg.norm(vector)
    if vNorm == 0.0:
        return vector
    else:
        return vector / vNorm

def calculate_angle_in_corner(p0, p1, p2):
    p10 = normalize(p0 - p1)
    p12 = normalize(p2 - p1)

    dp = max(-1.0, min(1.0, np.dot(p10, p12)))

    return np.arccos(dp)

def calculate_angle_between_vectors(v0, v1):
    dp = max(-1.0, min(1.0, np.dot(v0, v1)))

    return np.arccos(dp)

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

