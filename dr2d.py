import math
import numpy as np
import pyrender
import triangle
import trimesh

VIEWPORT_SIZE=(1280,720)
SMALLEST_ANGLE = math.radians(30)
SIZE_THRESHOLD = 0.5

def generatePSLGv1():
    vertices = [
        # boundary
        (-1.0, -1.0),
        ( 1.0, -1.0),
        ( 1.5,  0.0),
        ( 1.0,  1.0),
        (-1.0,  1.0),
        # hole
        (-0.3, -0.3),
        ( 0.3, -0.3),
        ( 0.3,  0.3),
        (-0.3,  0.3),
    ]
    segments = [
        # boundary
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 0),
        # hole
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 5)
    ]
    return vertices, segments

def generatePSLGv2():
    vertices = [
        ( 0.0,  0.0),
        (-1.0,  0.0),
        (-2.0,  1.0),
        (-2.0, -1.5),
        (-1.0, -0.5),
        ( 1.0, -0.5),
        ( 2.0, -1.5),
        ( 2.0,  1.0),
        ( 1.0,  0.0)
    ]
    segments = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 0)
    ]
    return vertices, segments

def generatePSLGv3():
    vertices = [(0.0, 0.0)] # circle center
    segments = []
    # gen circle section polygon with specified degree
    max_angle = 270
    degree = 300
    radius = 5
    for i in range(0, degree):
        p = 2*np.pi * i/degree
        if p > np.deg2rad(max_angle):
            break
        if (2*np.pi * (i+1)/degree) > np.deg2rad(max_angle):
            p = np.deg2rad(max_angle)
        x = radius * np.cos(p)
        y = radius * np.sin(p)
        vertices.append((x, y))
    for i in range(0, len(vertices)):
        s0 = i
        s1 = (i+1) % len(vertices)
        segments.append((s0, s1))
    return vertices, segments

def triangulate(pslg):
    vertices, segments = pslg

    t = triangle.Triangle()

    pointsBoundary = vertices
    markersBoundary = [1] * len(pointsBoundary)
    pointsInner = []
    markersInner = [0] * len(pointsInner)

    points = pointsBoundary + pointsInner
    markers = markersBoundary + markersInner
    t.set_points(points, markers=markers)
    t.set_segments(segments)

    t.triangulate(mode='pzQ')

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

def calcAngle(a, b, c):
    BA = np.array(a) - np.array(b)
    BC = np.array(c) - np.array(b)
    dp = np.dot(normalize(BA), normalize(BC))
    dp = min(1.0, max(-1.0, dp))
    return np.arccos(dp)

def shapeTest(cdt, iTriangle):
    vertices, _, triangles = cdt
    t0, t1, t2 = triangles[iTriangle]
    v0 = vertices[t0]
    v1 = vertices[t1]
    v2 = vertices[t2]
    alpha = calcAngle(v0, v1, v2)
    beta = calcAngle(v1, v2, v0)
    gamma = calcAngle(v2, v0, v1)
    return min(alpha, beta, gamma) > SMALLEST_ANGLE

def calcSize(a, b, c):
    AB = np.array(b) - np.array(a)
    AC = np.array(c) - np.array(a)
    size = np.cross(AB, AC) / 2.0
    assert size >= 0
    return size

def sizeTestv1(cdt, iTriangle):
    vertices, _, triangles = cdt
    t0, t1, t2 = triangles[iTriangle]
    v0 = vertices[t0]
    v1 = vertices[t1]
    v2 = vertices[t2]
    size = calcSize(v0, v1, v2)
    return size <= SIZE_THRESHOLD, size

def sizeTestv2(cdt, iTriangle):
    vertices, _, triangles = cdt
    t0, t1, t2 = triangles[iTriangle]
    v0 = vertices[t0]
    v1 = vertices[t1]
    v2 = vertices[t2]
    size = calcSize(v0, v1, v2)
    centerOfGravity = (np.array(v0) + np.array(v1) + np.array(v2)) / 3.0
    lowerBound = 0.001
    cogDist = max(lowerBound, np.linalg.norm(centerOfGravity))
    return 2*size <= (cogDist**2), size

def sizeTest(cdt, iTriangle):
    return sizeTestv1(cdt, iTriangle)

def findLargestFailingTriangle(cdt):
    vertices, segments, triangles = cdt
    delta = -1
    deltaSize = 0.0
    for iTriangle in range(0, len(triangles)):
        (t0, t1, t2) = triangles[iTriangle]
        wellShaped = shapeTest(cdt, iTriangle)
        wellSized, size = sizeTest(cdt, iTriangle)
        if wellShaped and wellSized:
            pass
        else:
            if size > deltaSize:
                deltaSize = size
                delta = iTriangle
    return delta

def calcCircumcenter(cdt, delta):
    vertices, _, triangles = cdt
    # https://stackoverflow.com/questions/56224824/how-do-i-find-the-circumcenter-of-the-triangle-using-python-without-external-lib
    t0, t1, t2 = triangles[delta]
    ax, ay = vertices[t0]
    bx, by = vertices[t1]
    cx, cy = vertices[t2]
    d = 2 * (ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
    ux = ((ax*ax + ay*ay) * (by-cy) + (bx*bx + by*by) * (cy-ay) + (cx*cx + cy*cy) * (ay-by)) / d
    uy = ((ax*ax + ay*ay) * (cx-bx) + (bx*bx + by*by) * (ax-cx) + (cx*cx + cy*cy) * (bx-ax)) / d
    return (ux, uy)

# https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
def onSegment(p, q, r):
    p_x, p_y = p
    q_x, q_y = q
    r_x, r_y = r
    return q_x <= max(p_x, r_x) and q_x >= min(p_x, r_x) and \
           q_y <= max(p_y, r_y) and q_y >= min(p_y, r_y)
def orientation(p, q, r):
    p_x, p_y = p
    q_x, q_y = q
    r_x, r_y = r
    val = (q_y - p_y) * (r_x - q_x) - (q_x - p_x) * (r_y - q_y)
    if val == 0:
        return 0
    else:
        return 1 if val > 0 else 2
def segmentsIntersect(seg0, seg1):
    p1, q1 = seg0
    p2, q2 = seg1
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)
    if o1 != o2 and o3 != o4:
        return True
    if o1 == 0 and onSegment(p1, p2, q1):
        return True
    if o2 == 0 and onSegment(p1, q2, q1):
        return True
    if o3 == 0 and onSegment(p2, p1, q2):
        return True
    if o4 == 0 and onSegment(p2, q1, q2):
        return True
    return False

def pointsAllClose(p0, p1):
    return np.allclose(np.array(p0), np.array(p1))

def oneSegmentEndpointAllClose(seg, p):
    return pointsAllClose(seg[0], p) or pointsAllClose(seg[1], p)

def getIntersectingSegmentIdHelper(cdt, v, c):
    vertices, segments, triangles = cdt
    for iSegment in range(0, len(segments)):
        (s0, s1) = segments[iSegment]
        sv0_sv1 = (vertices[s0], vertices[s1])
        # check for common endpoint
        assert not oneSegmentEndpointAllClose(sv0_sv1, c)
        if oneSegmentEndpointAllClose(sv0_sv1, v):
            continue
        # check for intersection
        if segmentsIntersect(sv0_sv1, (v, c)):
            return iSegment
        else:
            continue
    return -1

def getIntersectingSegmentId(cdt, delta, c):
    vertices, segments, triangles = cdt
    t0, t1, t2 = triangles[delta]
    v0 = vertices[t0]
    v1 = vertices[t1]
    v2 = vertices[t2]
    hit = getIntersectingSegmentIdHelper(cdt, v0, c)
    if (hit >= 0):
        return hit
    hit = getIntersectingSegmentIdHelper(cdt, v1, c)
    if (hit >= 0):
        return hit
    hit = getIntersectingSegmentIdHelper(cdt, v2, c)
    return hit

def calcHalfwayPoint(v0, v1):
    v0_x, v0_y = v0
    v1_x, v1_y = v1
    vh_x = v0_x + ((v1_x - v0_x) / 2.0)
    vh_y = v0_y + ((v1_y - v0_y) / 2.0)
    return (vh_x, vh_y)

def isSegmentVertex(cdt, iVertex):
    vertices, segments, triangles = cdt
    for s in segments:
        s0, s1 = s
        if s0 == iVertex or s1 == iVertex:
            return True
        else:
            continue
    return False

def findEncroachingInnerVertexId(cdt, vh, enchroachDist):
    vertices, segments, triangles = cdt
    for iVertex in range(0, len(vertices)):
        if isSegmentVertex(cdt, iVertex):
            continue
        v = vertices[iVertex]
        if v == vh: # avoid self-encroaching
            continue
        dist = np.linalg.norm(np.array(v)-np.array(vh))
        if dist < enchroachDist:
            return iVertex
        else:
            pass
    return -1

def removeVertex(cdt, removeId):
    vertices, segments, triangles = cdt
    del vertices[removeId]
    triangles.clear()
    for iSegment in range(0, len(segments)):
        s0, s1 = segments[iSegment]
        if s0 > removeId:
            s0 -= 1
        if s1 > removeId:
            s1 -= 1
        segments[iSegment] = s0, s1

def splitSegment(cdt, iSegment):
    vertices, segments, triangles = cdt
    (s0, s1) = segments[iSegment]
    # calculate and instert new halfway point
    v0 = vertices[s0]
    v1 = vertices[s1]
    vh = calcHalfwayPoint(v0, v1)
    sh = len(vertices)
    vertices.append(vh)
    # replace (s0, s1) segment by (s0, sh) and (sh , s1) segments
    del segments[iSegment]
    segments.append((s0, sh))
    segments.append((sh, s1))
    # remove all encroaching vertices
    enchroachDist = np.linalg.norm(np.array(vh)-np.array(v0))
    removeId =  findEncroachingInnerVertexId(cdt, vh, enchroachDist)
    while removeId >= 0:
        removeVertex(cdt, removeId)
        removeId =  findEncroachingInnerVertexId(cdt, vh, enchroachDist)
    return vertices, segments, triangles

def chew93(pslg):
    vertices, segments = pslg

    # step 1: build initial CDT from PSLG
    cdt = vertices, segments, triangulate(pslg)

    # step 2+3: find largest triangle that fails shape ans size criteria
    delta = findLargestFailingTriangle(cdt)

    iteration_counter = 0
    while delta >= 0:
        iteration_counter += 1
        if iteration_counter > 1 and (iteration_counter-1) % 100 == 0:
            print(iteration_counter-1, 'iterations processed...')
        # render(cdt)#TODO remove
        # step 4: travel from any triangle vertex to c and return hit segment id
        c = calcCircumcenter(cdt, delta)
        hit = getIntersectingSegmentId(cdt, delta, c)
        if hit < 0: # step 5: no segment was hit, insert c
            vertices.append(c)
        else: # step 6: split hit segment
            cdt = splitSegment(cdt, hit)
        # update cdt and delta
        cdt = vertices, segments, triangulate((vertices, segments))
        delta = findLargestFailingTriangle(cdt)
    print('converged after', iteration_counter, 'iterations')

    return cdt

def trimeshFromCdt(cdt):
    vertices = [list(v)+[-1.0] for v in cdt[0]]
    triangles = [list(t) for t in cdt[2]]
    mesh = trimesh.Trimesh(vertices=vertices, faces=triangles)

    return mesh

def render(cdt):
    tm = trimeshFromCdt(cdt)
    mesh = pyrender.Mesh.from_trimesh(tm, smooth=False)
    cam = pyrender.OrthographicCamera(xmag=5.0, ymag=5.0)
    pos = np.eye(4)
    scene = pyrender.Scene()
    scene.add(mesh)
    scene.add(cam, pose=pos)
    pyrender.Viewer(scene, VIEWPORT_SIZE, all_wireframe=True, cull_faces=False, use_perspective_cam=False)

if __name__ == '__main__':
    pslg = generatePSLGv3()
    cdt = chew93(pslg)
    render(cdt)
