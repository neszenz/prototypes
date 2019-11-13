""" Mesh1D class
This class represents the results of the sampling process.

name: path to the sampled step file
vertices: a complete list of all samples as SuperVertices
face_meshes: list of indexed face meshes
                        "=> tuple of...
                           ...a index loop of the outer boundary
                           ...a list of index loops for all inner boundaries
                           ...surface type string
bounding_box2D: for uv coordinates and -face_id as z
bounding_box3D: for xyz coordinates
"""
import numpy as np

BOUNDING_BOX_DEFAULT = (
    float('inf'), float('-inf'),
    float('inf'), float('-inf'),
    float('inf'), float('-inf')
)

FILE_EXTENSION = '.mesh1D'

class SuperVertex:
    def __init__(self, x=0.0, y=0.0, z=0.0, u=0.0, v=0.0, face_id=0):
        self.x = x
        self.y = y
        self.z = z
        self.u = u
        self.v = v
        self.face_id = face_id

    def UV_vec2(self):
        return np.array([self.u, self.v])
    def UV_vec3(self):
        return np.array([self.u, self.v, -(self.face_id)])
    def XYZ_vec3(self):
        return np.array([self.x, self.y, self.z])

    def allclose_UV(self, other):
        return np.allclose(self.UV_vec2(), other.UV_vec2())
    def allclose_XYZ(self, other):
        return np.allclose(self.XYZ_vec3(), other.XYZ_vec3())
    def __eq__(self, other):
        return self.allclose_UV(other) and self.allclose_XYZ(other) and self.face_id == other.face_id

    def __str__(self):
        return '(' + str(self.XYZ_vec3()) + ', ' + str(self.UV_vec2()) + ', ' + str(self.face_id) + ')'
    def __repr__(self):
        return self.__str__()

class Mesh1D:
    def __init__(self, name, vertices, face_meshes):
        self.name = name
        self.vertices = vertices
        assert all(isinstance(v, SuperVertex) for v in self.vertices)
        self.face_meshes = face_meshes
        self.bounding_box2D = BOUNDING_BOX_DEFAULT
        self.bounding_box3D = BOUNDING_BOX_DEFAULT
        self.update_bb2D()
        self.update_bb3D()

    def number_of_faces(self):
        return len(self.face_meshes)

    def update_bb2D(self):
        x_min, x_max, y_min, y_max, z_min, z_max = BOUNDING_BOX_DEFAULT
        for v in self.vertices:
            if v.u < x_min:
                x_min = v.u
            if v.u > x_max:
                x_max = v.u
            if v.v < y_min:
                y_min = v.v
            if v.v > y_max:
                y_max = v.v
            if -(v.face_id) < z_min:
                z_min = -(v.face_id)
            if -(v.face_id) > z_max:
                z_max = -(v.face_id)
        self.bounding_box2D = x_min, x_max, y_min, y_max, z_min, z_max
    def get_bb2D_size(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box2D
        return (x_max-x_min, y_max-y_min, z_max-z_min)
    def get_bb2D_size_factor(self):
        return np.linalg.norm(self.get_bb2D_size())
    def get_bb2D_center(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box2D
        cog = np.array((x_max-x_min, y_max-y_min, z_max-z_min)) / 2
        return np.array((x_min, y_min, z_min)) + cog

    def update_bb3D(self):
        x_min, x_max, y_min, y_max, z_min, z_max = BOUNDING_BOX_DEFAULT
        for v in self.vertices:
            if v.x < x_min:
                x_min = v.x
            if v.x > x_max:
                x_max = v.x
            if v.y < y_min:
                y_min = v.y
            if v.y > y_max:
                y_max = v.y
            if v.z < z_min:
                z_min = v.z
            if v.z > z_max:
                z_max = v.z
        self.bounding_box3D = x_min, x_max, y_min, y_max, z_min, z_max
    def get_bb3D_size(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box3D
        return (x_max-x_min, y_max-y_min, z_max-z_min)
    def get_bb3D_size_factor(self):
        return np.linalg.norm(self.get_bb3D_size())
    def get_bb3D_center(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box3D
        cog = np.array((x_max-x_min, y_max-y_min, z_max-z_min)) / 2
        return np.array((x_min, y_min, z_min)) + cog

    def get_face_type(self, index):
        _, _, face_type = self.face_meshes[index]
        return face_type
