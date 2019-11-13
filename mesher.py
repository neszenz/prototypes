""" Chew93_Surface mesher for STEP files using sampler.py w/ OpenCASCADE
Contains functionality to generate a 2-dimensional mesh from a given path to a
STEP file.

To gain a foundation for the meshing process, sampler.py is used to generate a
1-dimensional mesh from the STEP file which is than iterated upon. Here, the
first step is to compute constrained Delaunay Triangulations w/ pytriangle and
then apply my implementation of Chew's Surface Delaunay Refinement algorithm.
"""
import numpy as np
import triangle as shewchuck_triangle

import paths
import sampler

import pyrender #TODO remove
import trimesh #TODO remove
VIEWPORT_SIZE=(1280,720)
ORTHO_CAM = [VIEWPORT_SIZE[0]/10, VIEWPORT_SIZE[1]/10]

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_PATH = paths.PATH_98613
sampler.NUMBER_OF_SAMPLES = 20

def trimeshFromMesh(mesh): #TODO remove
    vertices, triangles = mesh
    vertices3D = [[sv.x, sv.y, sv.z] for sv in vertices]
    triangles = [list(t) for t in mesh[1]]
    mesh = trimesh.Trimesh(vertices=vertices3D, faces=triangles)
    return mesh

def renderMesh(mesh): #TODO remove
    tm = trimeshFromMesh(mesh)
    sceneMesh = pyrender.Mesh.from_trimesh(tm, smooth=False)
    scene = pyrender.Scene()
    scene.add(sceneMesh)
    pyrender.Viewer(scene, VIEWPORT_SIZE, all_wireframe=True, cull_faces=False, use_perspective_cam=True)

def triangulate(mesh1D):
    def collect_vertex_indices_of_face(face_mesh):
        vertices = []
        outer_wire, inner_wires, _ = face_mesh
        vertices += outer_wire
        for inner_wire in inner_wires:
            vertices += inner_wire
        return vertices
    def vertices_from_vertex_indices(mesh1D, vertex_indices):
        vertices = []
        for vertex_index in vertex_indices:
            vertices.append(mesh1D.vertices[vertex_index])
        return vertices
    def triangulate_face(vertices, vertex_offset):
        triangles = []
        t = shewchuck_triangle.Triangle()
        pointsBoundary = [[sv.u, sv.v] for sv in vertices]
        markersBoundary = [1] * len(pointsBoundary)
        pointsInner = []
        markersInner = [0] * len(pointsInner)
        points = pointsBoundary + pointsInner
        markers = markersBoundary + markersInner
        t.set_points(points, markers=markers)
        t.triangulate(mode='zQ')
        for triNode in t.get_triangles():
            ([i0, i1, i2], neighbors, attri) = triNode
            i0 += vertex_offset
            i1 += vertex_offset
            i2 += vertex_offset
            triangles.append((i0, i1, i2))
        return triangles
    vertices = []
    triangles = []
    for face_mesh in mesh1D.face_meshes:
        face_vertex_indices = collect_vertex_indices_of_face(face_mesh)
        face_vertices = vertices_from_vertex_indices(mesh1D, face_vertex_indices)
        triangles+= triangulate_face(face_vertices, len(vertices))
        vertices += face_vertices
    return vertices, triangles

if __name__ == '__main__':
    simple_sampler = sampler.factory(sampler.SAMPLER_TYPE.SIMPLE)
    mesh1D = simple_sampler(INPUT_PATH)
    sampler.write_mesh_to_file(mesh1D)
    mesh2D = triangulate(mesh1D)
    renderMesh(mesh2D) #TODO remove
