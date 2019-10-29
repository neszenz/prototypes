import math
import numpy as np
import pyrender
import random
import triangle
import trimesh

VIEWPORT_SIZE=(1280,720)
DISTRIBUTION_WINDOW = [VIEWPORT_SIZE[0]/10, VIEWPORT_SIZE[1]/10]
ORTHO_CAM = [VIEWPORT_SIZE[0]/10, VIEWPORT_SIZE[1]/10]
NUM_SAMPLES = 5

def generatePointSet():
    vertices = []
    for i in range(0, NUM_SAMPLES):
        x = (2.0 * DISTRIBUTION_WINDOW[0] * random.random()) - DISTRIBUTION_WINDOW[0]
        y = (2.0 * DISTRIBUTION_WINDOW[1] * random.random()) - DISTRIBUTION_WINDOW[1]
        vertices.append([x, y])
    return vertices

def generatePointSetv2():
    x_norm = DISTRIBUTION_WINDOW[0]
    y_norm = DISTRIBUTION_WINDOW[1]
    vertices = [
        [-x_norm,-y_norm],
        [ x_norm,-y_norm],
        [ x_norm, y_norm],
        [-x_norm, y_norm],
        [0,0]
    ]
    return vertices

def triangulate(vertices):
    t = triangle.Triangle()

    pointsBoundary = vertices
    markersBoundary = [1] * len(pointsBoundary)
    pointsInner = []
    markersInner = [0] * len(pointsInner)

    points = pointsBoundary + pointsInner
    markers = markersBoundary + markersInner
    t.set_points(points, markers=markers)

    t.triangulate(mode='zQ')

    triangles = []
    for triNode in t.get_triangles():
        ([i0, i1, i2], neighbors, attri) = triNode
        triangles.append((i0, i1, i2))

    return triangles

def normalize(vector):
    vNorm = np.linalg.norm(vector)
    if vNorm == 0.0:
        return vector
    else:
        return vector / vNorm

def calcCircumcenter(mesh, delta):
    vertices, triangles = mesh
    # https://stackoverflow.com/questions/56224824/how-do-i-find-the-circumcenter-of-the-triangle-using-python-without-external-lib
    t0, t1, t2 = triangles[delta]
    ax, ay = vertices[t0]
    bx, by = vertices[t1]
    cx, cy = vertices[t2]
    d = 2 * (ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
    ux = ((ax*ax + ay*ay) * (by-cy) + (bx*bx + by*by) * (cy-ay) + (cx*cx + cy*cy) * (ay-by)) / d
    uy = ((ax*ax + ay*ay) * (cx-bx) + (bx*bx + by*by) * (ax-cx) + (cx*cx + cy*cy) * (bx-ax)) / d
    return (ux, uy)

def circumcentersOfMesh(mesh):
    vertices, triangles = mesh
    circumcenters = []
    for iTriangle in range(0, len(triangles)):
        cc = calcCircumcenter(mesh, iTriangle)
        circumcenters.append(cc)
    return circumcenters

def pointsAllClose(p0, p1):
    return np.allclose(np.array(p0), np.array(p1))

def trimeshFromMesh(mesh):
    vertices, triangles = mesh
    vertices3D = [list(v)+[-1.0] for v in vertices]
    triangles = [list(t) for t in mesh[1]]
    mesh = trimesh.Trimesh(vertices=vertices3D, faces=triangles)

    return mesh

def renderMesh(mesh):
    tm = trimeshFromMesh(mesh)
    sceneMesh = pyrender.Mesh.from_trimesh(tm, smooth=False)
    cam = pyrender.OrthographicCamera(xmag=ORTHO_CAM[0], ymag=ORTHO_CAM[1])
    pos = np.eye(4)
    scene = pyrender.Scene()
    scene.add(sceneMesh)
    scene.add(cam, pose=pos)
    pyrender.Viewer(scene, VIEWPORT_SIZE, all_wireframe=True, cull_faces=False, use_perspective_cam=False)

def renderMeshWithCircumcenters(mesh, circumcenters):
    tm = trimeshFromMesh(mesh)
    sceneMesh = pyrender.Mesh.from_trimesh(tm, smooth=False)
    scenePoints = pyrender.Mesh.from_points([[v[0], v[1], -1.0] for v in circumcenters])
    cam = pyrender.OrthographicCamera(xmag=ORTHO_CAM[0], ymag=ORTHO_CAM[1])
    pos = np.eye(4)
    scene = pyrender.Scene()
    scene.add(sceneMesh)
    scene.add(scenePoints)
    scene.add(cam, pose=pos)
    pyrender.Viewer(scene, VIEWPORT_SIZE, all_wireframe=True, cull_faces=False, use_perspective_cam=False, point_size=5)

def renderVertices(vertices):
    mesh = pyrender.Mesh.from_points([[v[0], v[1], -1.0] for v in vertices])
    cam = pyrender.OrthographicCamera(xmag=ORTHO_CAM[0], ymag=ORTHO_CAM[1])
    pos = np.eye(4)
    scene = pyrender.Scene()
    scene.add(mesh)
    scene.add(cam, pose=pos)
    pyrender.Viewer(scene, VIEWPORT_SIZE, all_wireframe=True, cull_faces=False, use_perspective_cam=False)

if __name__ == '__main__':
    vertices = generatePointSet()
    mesh = vertices, triangulate(vertices)
    circumcenters = circumcentersOfMesh(mesh)
    renderMeshWithCircumcenters(mesh, circumcenters)
