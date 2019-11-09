""" Mesh1D class
This class represents the results of the sampling process.

name: path to the sampled step file
vertices: a complete list of all sampled vertices
face_meshes: list of indexed face meshes
                        "=> tuple of...
                           ...a index loop of the outer boundary
                           ...a list of index loops for all inner boundaries
                           ...surface type string
bounding_box: tuple of coordinate value ranges computed based on vertices
"""
import numpy as np

BOUNDING_BOX_DEFAULT = (
    float('inf'), float('-inf'),
    float('inf'), float('-inf'),
    float('inf'), float('-inf')
)

FILE_EXTENSION = '.mesh1D'

class Mesh1D:
    def __init__(self, name, vertices, face_meshes):
        self.name = name
        self.vertices = vertices
        self.face_meshes = face_meshes
        self.bounding_box = BOUNDING_BOX_DEFAULT
        self.update_bb()

    def number_of_faces(self):
        return len(self.face_meshes)

    def update_bb(self):
        x_min, x_max, y_min, y_max, z_min, z_max = BOUNDING_BOX_DEFAULT
        for v in self.vertices:
            if v[0] < x_min:
                x_min = v[0]
            if v[0] > x_max:
                x_max = v[0]
            if v[1] < y_min:
                y_min = v[1]
            if v[1] > y_max:
                y_max = v[1]
            if v[2] < z_min:
                z_min = v[2]
            if v[2] > z_max:
                z_max = v[2]
        self.bounding_box = x_min, x_max, y_min, y_max, z_min, z_max

    def get_bb_size(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box
        return (x_max-x_min, y_max-y_min, z_max-z_min)

    def get_bb_size_factor(self):
        return np.linalg.norm(self.get_bb_size())

    def get_bb_center(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box
        cog = np.array((x_max-x_min, y_max-y_min, z_max-z_min)) / 2
        return np.array((x_min, y_min, z_min)) + cog

    def get_face_type(self, index):
        _, _, face_type = self.face_meshes[index]
        return face_type
