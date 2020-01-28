""" MeshkD class
This class is used to store the results of the sampling and meshing process.

name: e.g. source file path
face_meshes: list of face meshes
                        "=> tuple of...
                            ...list of super vertices
                            ...list of wire meshes (first is outer; after that inner)
                                        "=> loop-list of vertex indices
                            ...list of triangles
                                        "=> 3-tuple of vertex indices
                            ...TopoDS_FACE
bounding_box2D: for uv coordinates and negative face_id as w
bounding_box3D: for xyz coordinates
"""
import datetime
import os
import pickle
import numpy as np

from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCC.Core.GCPnts import GCPnts_AbscissaPoint
from OCC.Core.GeomAbs import GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_BezierSurface, GeomAbs_BSplineSurface, GeomAbs_SurfaceOfRevolution, GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface, GeomAbs_OtherSurface
from OCC.Core.gp import gp_Pnt
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_FORWARD, TopAbs_REVERSED

from util import *

BOUNDING_BOX_DEFAULT = (
    float('inf'), float('-inf'),
    float('inf'), float('-inf'),
    float('inf'), float('-inf')
)

SURFACE_TYPE_STRINGS = {
    GeomAbs_Plane :               'Plane',
    GeomAbs_Cylinder :            'Cylinder',
    GeomAbs_Cone :                'Cone',
    GeomAbs_Sphere :              'Sphere',
    GeomAbs_Torus :               'Torus',
    GeomAbs_BezierSurface :       'BezierSurface',
    GeomAbs_BSplineSurface :      'BSplineSurface',
    GeomAbs_SurfaceOfRevolution : 'SurfaceOfRevolution',
    GeomAbs_SurfaceOfExtrusion :  'SurfaceOfExtrusion',
    GeomAbs_OffsetSurface :       'OffsetSurface',
    GeomAbs_OtherSurface :        'OtherSurface'
}

## SuperVertex + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
class SuperVertex:
    def __init__(self, x=0.0, y=0.0, z=0.0, u=0.0, v=0.0):
        self.x = x
        self.y = y
        self.z = z
        self.u = u
        self.v = v

        self.face_id = 0
        self.face = None

        # None: inner vertex
        # len-1-list: vertex on boundary
        # len-2-list: vertex on boundary endpoint
        self.edges_with_p = None # tuple of TopoDS_EDGE and sample parameter

    def set_same_face_as(self, other):
        self.face = other.face
        self.face_id = other.face_id
        return

    def compute_halfway_on_shared_edge(sv0, sv1):
        def get_shared_edge(sv0, sv1):
            for sv0_edge_with_p in sv0.edges_with_p:
                sv0_edge, sv0_p = sv0_edge_with_p

                for sv1_edge_with_p in sv1.edges_with_p:
                    sv1_edge, sv1_p = sv1_edge_with_p

                    if sv0_edge == sv1_edge:
                        u0 = min(sv0_p, sv1_p)
                        u1 = max(sv0_p, sv1_p)

                        return sv0_edge, u0, u1

            return None, -1, -1
        assert sv0.face_id== sv1.face_id
        assert not sv0.edges_with_p is None
        assert not sv1.edges_with_p is None

        shared_edge, u0, u1 = get_shared_edge(sv0, sv1)
        if shared_edge is None:
            raise Exception('compute_halfway_on_shared_edge() error - no shared edge')
        curve = BRepAdaptor_Curve2d(shared_edge, sv0.face)
        u01_length = GCPnts_AbscissaPoint.Length(curve, u0, u1)

        abscissa_point = GCPnts_AbscissaPoint(curve, u01_length/2, u0)
        assert abscissa_point.IsDone()

        p = abscissa_point.Parameter()
        uv = curve.Value(p)
        xyz = BRepAdaptor_Surface(sv0.face).Value(uv.X(), uv.Y())

        sv_halfway = SuperVertex(x=xyz.X(), y=xyz.Y(), z=xyz.Z(), u=uv.X(), v=uv.Y())
        sv_halfway.face_id = sv0.face_id
        sv_halfway.face = sv0.face
        sv_halfway.edges_with_p = [(shared_edge, p)]
        if sv0.face.Orientation() == TopAbs_REVERSED:
            sv_halfway.reverse_u()

        return sv_halfway
    def compute_halfway(sv0, sv1):
        p0 = sv0.UV_vec2()
        p1 = sv1.UV_vec2()
        phw = p0 + ((p1-p0)/2)
        sv_halfway = SuperVertex(u=phw[0], v=phw[1])
        sv_halfway.face_id = sv0.face_id
        sv_halfway.face = sv0.face
        sv_halfway.project_to_XYZ()

        return sv_halfway

    def UV_vec2(self):
        return np.array([self.u, self.v])
    def UV_vec3(self):
        return np.array([self.u, self.v, -self.face_id])
    def XYZ_vec3(self):
        return np.array([self.x, self.y, self.z])

    def reverse_u(self):
        self.u = reverse_u(self.u, self.face)

    def project_to_UV(self):
        assert self.face != None
        surface = BRep_Tool.Surface(self.face)
        analysis_surface = ShapeAnalysis_Surface(surface)
        xyz = gp_Pnt(self.x, self.y, self.z)
        uv = analysis_surface.ValueOfUV(xyz, 0.0001)
        self.u = uv.X()
        self.v = uv.Y()
        if self.face.Orientation() == TopAbs_REVERSED:
            self.reverse_u()
    def project_to_XYZ(self):
        assert self.face != None
        surface = BRepAdaptor_Surface(self.face)
        if self.face.Orientation() == TopAbs_REVERSED:
            u_tmp = self.u
            self.reverse_u()
            xyz = surface.Value(self.u, self.v)
            self.u = u_tmp
        else:
            xyz = surface.Value(self.u, self.v)
        self.x = xyz.X()
        self.y = xyz.Y()
        self.z = xyz.Z()

    def allclose_UV(self, other):
        return np.allclose(self.UV_vec2(), other.UV_vec2())
    def allclose_XYZ(self, other):
        return np.allclose(self.XYZ_vec3(), other.XYZ_vec3())
    def __eq__(self, other):
        return self.allclose_UV(other) and self.allclose_XYZ(other) and self.face_id == other.face_id

    def __str__(self):
        return '(' + str(self.XYZ_vec3()) + ', ' + str(self.UV_vec2()) + ', ' + str(self.face_id) + ', ' + str(self.edges_with_p) + ')'
    def __repr__(self):
        return self.__str__()

## MeshkD  + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
class MeshkD:
    FILE_EXTENSION = '.meshkD'
    NV_INIT = 'initial vertices: '
    NV_REFI = 'refinement vertices: '
    NV_INST = 'inner vertex inserts: '
    NV_SPLT = 'segment splits: '
    NV_DELT = 'encroaching vertices deleted: '
    NT_INVK = 'travel test invokations: '
    NT_LOOP = 'total travel test loops: '
    NT_SHAD = 'ray shadow tests: '

    def __init__(self, name, face_meshes):
        self.name = name
        self.face_meshes = face_meshes
        self.bounding_box3D = BOUNDING_BOX_DEFAULT
        self.bounding_box2D = BOUNDING_BOX_DEFAULT
        self.reset_bounding_boxes()

        self.meta_blocks = []
        for face_mesh in self.face_meshes:
            vertices, _, _, _ = face_mesh
            meta_block = {
                MeshkD.NV_INIT:len(vertices),
                MeshkD.NV_REFI:0,
                MeshkD.NV_INST:0,
                MeshkD.NV_SPLT:0,
                MeshkD.NV_DELT:0,
                MeshkD.NT_INVK:0,
                MeshkD.NT_LOOP:0,
                MeshkD.NT_SHAD:0
            }
            self.meta_blocks.append(meta_block)

    def number_of_faces(self):
        return len(self.face_meshes)
    def get_face_type(self, index):
        _, _, _, face = self.face_meshes[index]
        return SURFACE_TYPE_STRINGS[BRepAdaptor_Surface(face).GetType()]

    def reset_bounding_boxes(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box3D
        u_min, u_max, v_min, v_max, w_min, w_max = self.bounding_box2D

        for face_mesh in self.face_meshes:
            vertices, _, _, _ = face_mesh

            for sv in vertices:
                x_min = sv.x if sv.x < x_min else x_min
                x_max = sv.x if sv.x > x_max else x_max
                y_min = sv.y if sv.y < y_min else y_min
                y_max = sv.y if sv.y > y_max else y_max
                z_min = sv.z if sv.z < z_min else z_min
                z_max = sv.z if sv.z > z_max else z_max
                u_min = sv.u if sv.u < u_min else u_min
                u_max = sv.u if sv.u > u_max else u_max
                v_min = sv.v if sv.v < v_min else v_min
                v_max = sv.v if sv.v > v_max else v_max
                w_min = -sv.face_id if -sv.face_id < w_min else w_min
                w_max = -sv.face_id if -sv.face_id > w_max else w_max

        self.bounding_box3D = x_min, x_max, y_min, y_max, z_min, z_max
        self.bounding_box2D = u_min, u_max, v_min, v_max, w_min, w_max

        return

    def get_bb2D_size(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box2D
        return (x_max-x_min, y_max-y_min, z_max-z_min)
    def get_bb2D_size_factor(self):
        return np.linalg.norm(self.get_bb2D_size())
    def get_bb2D_center(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box2D
        cog = np.array((x_max-x_min, y_max-y_min, z_max-z_min)) / 2
        return np.array((x_min, y_min, z_min)) + cog

    def get_bb3D_size(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box3D
        return (x_max-x_min, y_max-y_min, z_max-z_min)
    def get_bb3D_size_factor(self):
        return np.linalg.norm(self.get_bb3D_size())
    def get_bb3D_center(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box3D
        cog = np.array((x_max-x_min, y_max-y_min, z_max-z_min)) / 2
        return np.array((x_min, y_min, z_min)) + cog

## I/O + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
def write_to_file(meshkD, output_dir):
    timestamp = datetime.datetime.now()
    name = timestamp.strftime('%y%m%d_%H%M%S') + MeshkD.FILE_EXTENSION
    path = os.path.join(output_dir, name)

    if os.path.exists(path):
        name = timestamp.strftime('%y%m%d_%H%M%S%f') + MeshkD.FILE_EXTENSION
        path = os.path.join(output_dir, name)
    assert not os.path.exists(path)

    print('>> write to file \"', path, '\"', sep='')
    os.makedirs(output_dir, exist_ok=True)
    pickle.dump(meshkD, open(path, 'wb'))

    return
def load_from_file(path):
    return pickle.load(open(path, 'rb'))
