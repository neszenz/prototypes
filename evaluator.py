from meshkD import SuperVertex, MeshkD, load_from_file
from util import *

# metric configuration
METRIC_ANGLES_HIST = False
METRIX_XI_HIST = False

# number of bins for metric histograms
METRIC_ANGLES_NOB = 53
METRIX_XI_NOB = 32

# __main__ config
INPUT_PATH = '../results/tmp/200118_204525.meshkD'

def basic_information(meshkD):
    print('\n -- basic information --')

    print('based on file ', end='')
    print(meshkD.name)

    num_of_faces = meshkD.number_of_faces()
    for i_face in range(meshkD.number_of_faces()):
        vertices, wire_meshes, triangles, face = meshkD.face_meshes[i_face]
        print('face_mesh ', i_face+1, ' of ', num_of_faces, ':', sep='')

        print('|   vertices: ', end='')
        print(len(vertices))

        num_of_wires = len(wire_meshes)
        for i_wire in range(len(wire_meshes)):
            wire_mesh = wire_meshes[i_wire]
            print('|   wire_mesh ', i_wire+1, ' of ', num_of_wires, ':', sep='')

            print('|   |   vertices: ', end='')
            print(len(wire_mesh))

        print('|   triangles: ', end='')
        print(len(triangles))

        print('|   surface type: ', end='')
        print(meshkD.get_face_type(i_face))

    return

def triangle_angle_metric(meshkD):
    def calculate_all_angles(vertices, triangles):
        angles = []
        min_angle = float('inf')
        max_angle = float('-inf')

        for (i0, i1, i2) in triangles:
            sv0 = vertices[i0]
            sv1 = vertices[i1]
            sv2 = vertices[i2]

            a0 = calculate_angle_in_corner(sv2.XYZ_vec3(), sv0.XYZ_vec3(), sv1.XYZ_vec3())
            a1 = calculate_angle_in_corner(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3())
            a2 = calculate_angle_in_corner(sv1.XYZ_vec3(), sv2.XYZ_vec3(), sv0.XYZ_vec3())

            min_a = min(a0, a1, a2)
            if min_a < min_angle:
                min_angle = min_a
            max_a = max(a0, a1, a2)
            if max_a > max_angle:
                max_angle = max_a

            angles += [a0, a1, a2]

        return angles, min_angle, max_angle
    print('\n -- trangle angles --')

    num_of_faces = meshkD.number_of_faces()
    for i_face in range(meshkD.number_of_faces()):
        vertices, _, triangles, _ = meshkD.face_meshes[i_face]
        print('face_mesh ', i_face+1, ' of ', num_of_faces, ':', sep='')

        angles, min_angle, max_angle = calculate_all_angles(vertices, triangles)
        print('|   range of angles: ', end='')
        print('[', round(np.rad2deg(min_angle), 2), end='')
        print(', ', round(np.rad2deg(max_angle), 2), ']', sep='')

        if METRIC_ANGLES_HIST:
            angles_degree = [np.rad2deg(a) for a in angles]
            plot_histogram(angles_degree, METRIC_ANGLES_NOB)

    return

def normalized_area_metric(meshkD):
    def calculate_xi_for_all_triangles(vertices, triangles):
        def calculate_xi_for_triangle(p0, p1, p2):
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
        xi_list = []
        min_xi = float('inf')
        max_xi = float('-inf')

        for (i0, i1, i2) in triangles:
            sv0 = vertices[i0]
            sv1 = vertices[i1]
            sv2 = vertices[i2]

            xi = calculate_xi_for_triangle(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3())
            if xi < min_xi:
                min_xi = xi
            if xi > max_xi:
                max_xi = xi

            xi_list.append(xi)

        return xi_list, min_xi, max_xi
    print('\n -- normalized area metric (xi-metric) --')

    num_of_faces = meshkD.number_of_faces()
    for i_face in range(meshkD.number_of_faces()):
        vertices, _, triangles, _ = meshkD.face_meshes[i_face]
        print('face_mesh ', i_face+1, ' of ', num_of_faces, ':', sep='')

        xi_list, min_xi, max_xi = calculate_xi_for_all_triangles(vertices, triangles)
        print('|   range of xi: ', end='')
        print('[', min_xi, ', ', max_xi, ']', sep='')

        if METRIX_XI_HIST:
            xi_lower_bound = 0.0
            xi_upper_bound = np.sqrt(3) / 4
            plot_histogram(xi_list, METRIX_XI_NOB, xi_lower_bound, xi_upper_bound)

def evaluate(path):
    print('>> evaluating file', path)
    meshkD = load_from_file(path)

    basic_information(meshkD)

    triangle_angle_metric(meshkD)

    normalized_area_metric(meshkD)

    print('')

if __name__ == '__main__':
    evaluate(INPUT_PATH)
