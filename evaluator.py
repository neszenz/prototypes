import os.path
import pathlib

from meshkD import SuperVertex, MeshkD, load_from_file
from util import *

# metric's constants, flags and configuration
N_DECIMAL_PLACES = 2
PLOT_SIZE = 10

METRIC_ANGLES     = 0
METRIC_MIN_ANGLES = 1
METRIC_XI         = 2
METRIC_RATIO      = 3

METRIC_NAMES = {
    METRIC_ANGLES :     'histogram of all angles',
    METRIC_MIN_ANGLES : 'histogram of minimum angles',
    METRIC_XI :         'histogram of xi metric',
    METRIC_RATIO :      'histogram of radius-edge ratio'
}
METRIC_HIST_FLAGS = {
    METRIC_ANGLES :     False,
    METRIC_MIN_ANGLES : False,
    METRIC_XI :         True,
    METRIC_RATIO :      False
}
METRIC_NOB = {
    METRIC_ANGLES :     70,
    METRIC_MIN_ANGLES : 70,
    METRIC_XI :         80,
    METRIC_RATIO :      80
}
METRIC_RANGES = {
    METRIC_ANGLES :     (0.0, np.rad2deg(np.pi)),
    METRIC_MIN_ANGLES : (0.0, np.rad2deg(np.pi)),
    METRIC_XI :         (0.0, np.sqrt(3) / 4),
    METRIC_RATIO :      (1.0/np.sqrt(3), float('inf'))
}

# __main__ config
INPUT_DIR = '../results/thesis_figures/element_quality/cylinder_progression/'
INPUT_PATH = '../results/tmp/200229_194956.meshkD'
INDENT = '|   '

def basic_information(mesh):
    def print_meta_block(meta_block, indentaion):
        print(INDENT + 'meta_block:')
        for item in meta_block:
            print(INDENT + INDENT + item + str(meta_block[item]))
    print('\n -- basic information --')

    print('based on file ', end='')
    print(mesh.name)

    num_of_faces = mesh.number_of_faces()
    for i_face in range(mesh.number_of_faces()):
        vertices, wire_meshes, triangles, face = mesh.face_meshes[i_face]
        print('face_mesh ', i_face+1, ' of ', num_of_faces, ':', sep='')

        print(INDENT + 'vertices: ', end='')
        print(len(vertices))

        print_meta_block(mesh.meta_blocks[i_face], INDENT)

        num_of_wires = len(wire_meshes)
        for i_wire in range(len(wire_meshes)):
            wire_mesh = wire_meshes[i_wire]
            print(INDENT + 'wire_mesh ', i_wire+1, ' of ', num_of_wires, ':', sep='')

            print(INDENT + INDENT + 'vertices: ', end='')
            print(len(wire_mesh))

        print(INDENT + 'triangles: ', end='')
        print(len(triangles))

        print(INDENT + 'surface type: ', end='')
        print(mesh.get_face_type(i_face))

    return

def calculate_metric_values(vertices, triangles, metric):
    def calculate_all_angles(vertices, triangles):
        angle_list = []
        min_angle = float('inf')
        max_angle = float('-inf')

        for (i0, i1, i2) in triangles:
            sv0 = vertices[i0]
            sv1 = vertices[i1]
            sv2 = vertices[i2]

            a0 = calculate_angle_in_corner(sv2.XYZ_vec3(), sv0.XYZ_vec3(), sv1.XYZ_vec3())
            a1 = calculate_angle_in_corner(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3())
            a2 = calculate_angle_in_corner(sv1.XYZ_vec3(), sv2.XYZ_vec3(), sv0.XYZ_vec3())

            a0 = np.rad2deg(a0)
            a1 = np.rad2deg(a1)
            a2 = np.rad2deg(a2)

            min_a = min(a0, a1, a2)
            max_a = max(a0, a1, a2)

            if min_a < min_angle:
                min_angle = min_a
            if max_a > max_angle:
                max_angle = max_a

            angle_list += [a0, a1, a2]

        return angle_list, min_angle, max_angle
    def calculate_min_angles(vertices, triangles):
        min_angle_list = []
        min_angle = float('inf')
        max_angle = float('-inf')

        for (i0, i1, i2) in triangles:
            sv0 = vertices[i0]
            sv1 = vertices[i1]
            sv2 = vertices[i2]

            a0 = calculate_angle_in_corner(sv2.XYZ_vec3(), sv0.XYZ_vec3(), sv1.XYZ_vec3())
            a1 = calculate_angle_in_corner(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3())
            a2 = calculate_angle_in_corner(sv1.XYZ_vec3(), sv2.XYZ_vec3(), sv0.XYZ_vec3())

            a0 = np.rad2deg(a0)
            a1 = np.rad2deg(a1)
            a2 = np.rad2deg(a2)

            min_a = min(a0, a1, a2)
            if min_a < min_angle:
                min_angle = min_a
            if min_a > max_angle:
                max_angle = min_a

            min_angle_list.append(min_a)

        return min_angle_list, min_angle, max_angle
    def calculate_xi_for_all_triangles(vertices, triangles):
        xi_list = []
        min_xi = float('inf')
        max_xi = float('-inf')

        for (i0, i1, i2) in triangles:
            sv0 = vertices[i0]
            sv1 = vertices[i1]
            sv2 = vertices[i2]

            xi = calculate_xi(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3())
            if xi < min_xi:
                min_xi = xi
            if xi > max_xi:
                max_xi = xi

            xi_list.append(xi)

        return xi_list, min_xi, max_xi
    def calculate_ratio_for_all_triangles(vertices, triangles):
        ratio_list = []
        min_ratio = float('inf')
        max_ratio = float('-inf')

        for (i0, i1, i2) in triangles:
            sv0 = vertices[i0]
            sv1 = vertices[i1]
            sv2 = vertices[i2]

            ratio = calculate_radius_edge_ratio_v2(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3())
            if ratio < min_ratio:
                min_ratio = ratio
            if ratio > max_ratio:
                max_ratio = ratio

            ratio_list.append(ratio)

        return ratio_list, min_ratio, max_ratio
    if metric == METRIC_ANGLES:
        return calculate_all_angles(vertices, triangles)
    elif metric == METRIC_MIN_ANGLES:
        return calculate_min_angles(vertices, triangles)
    elif metric == METRIC_XI:
        return calculate_xi_for_all_triangles(vertices, triangles)
    elif metric == METRIC_RATIO:
        return calculate_ratio_for_all_triangles(vertices, triangles)
    else:
        raise Exception('calculate_metric_values() error - unknown metric')

def evaluate_metric(mesh, metric):
    print('\n -- ' + METRIC_NAMES[metric] + ' --')

    num_of_faces = mesh.number_of_faces()
    for i_face in range(num_of_faces):
        vertices, _, triangles, _ = mesh.face_meshes[i_face]
        print('face_mesh ', i_face+1, ' of ', num_of_faces, ':', sep='')

        values, min_value, max_value = calculate_metric_values(vertices, triangles, metric)
        print(INDENT + 'range: ', end='')
        print('[', round(min_value, N_DECIMAL_PLACES), sep='', end='')
        print(', ', round(max_value, N_DECIMAL_PLACES), ']', sep='')

    if METRIC_HIST_FLAGS[metric]:
        file_name = str(metric)+'/'
        file_name += os.path.basename(INPUT_PATH)[:-len(MeshkD.FILE_EXTENSION)]
        iter_str = str(101) if cnt == 21 else str(cnt*5)
        title = METRIC_NAMES[metric]+' ('+iter_str+' iterations)'
        plot_histogram(values, title, METRIC_NOB[metric], METRIC_RANGES[metric], PLOT_SIZE, file_name)

    return

def evaluate(path):
    print('>> evaluating file', path)
    mesh = load_from_file(path)

    # basic_information(mesh)

    evaluate_metric(mesh, METRIC_ANGLES)
    evaluate_metric(mesh, METRIC_MIN_ANGLES)
    evaluate_metric(mesh, METRIC_XI)
    evaluate_metric(mesh, METRIC_RATIO)

    print('')
    return

if __name__ == '__main__':
    paths = pathlib.Path(INPUT_DIR).rglob('*'+MeshkD.FILE_EXTENSION)
    cnt = 0
    for path in sorted(paths):
        INPUT_PATH = path
        evaluate(INPUT_PATH)
        cnt += 1
