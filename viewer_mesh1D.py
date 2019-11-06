import os
import glob
import numpy as np
import pickle

import pyglet
from pyglet.gl import *
from pyglet.window import key
from pyglet.window import mouse

import mesh1D
from mesh1D import Mesh1D

## config constants  = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_DIR = 'tmp'
DEFAULT_X = 1280
DEFAULT_Y = 720
ZOOM_DEFAULT = 100
ZOOM_STEP = ZOOM_DEFAULT / 10
SAVE_ZONE_FACTOR = 0.8 # scale weight for model size to window ratio
ROTATION_DEFAULT = 0
ROTATION_STEP_SLOW = 1
ROTATION_STEP_FAST = 5
TILT_DEFAULT = 0

## global variables  = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
meshes = []
window = pyglet.window.Window(DEFAULT_X, DEFAULT_Y,
                              fullscreen=False,
                              resizable=True,
                              vsync=True)
window_caption = 'Delaunay Surface Mesh Refinement'

flags = {
    'draw_bounding_box' : True,
    'draw_vertices' : True,
    'draw_mesh' : True,
    'arrow_mode' : False, # draw lines as arrows; only in individual face mode
    'mesh_drawn_since_last_face_index_change' : False # avoids skipping faces while flipping through
}
draw_mesh_index = 0
draw_face_index = 0 # 0: all faces; >1: select individual face
direction_of_flight = np.array([0.0, 0.0, 1.0])
zoom = ZOOM_DEFAULT
rotation = ROTATION_DEFAULT
tilt = TILT_DEFAULT

def load_meshes_from_files():
    global meshes
    meshes.clear()
    files = sorted(glob.glob(os.path.join(INPUT_DIR, '**/*'+mesh1D.FILE_EXTENSION)), reverse=True)
    for path in files:
        loaded_mesh1D = pickle.load(open(path, 'rb'))
        loaded_mesh1D.name = path
        meshes.append(loaded_mesh1D)

def resetProjectionMatrix(width):
    w = width
    h = w / (window.width/window.height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(-w/2, w/2, -h/2, h/2, -w, w)
    glMatrixMode(GL_MODELVIEW)

@window.event
def on_resize(width, height):
    glViewport(0, 0, width, height)
    resetProjectionMatrix(width)
    return pyglet.event.EVENT_HANDLED

def apply_to_rotation(value):
    global rotation
    if rotation + value < 0:
        rotation = 360 + value
    elif rotation + value > 360:
        rotation = value
    elif rotation + value == 360:
            rotation = 0
    else:
        rotation += value

@window.event
def on_text_motion(motion):
    def apply_to_zoom(value):
        global zoom
        if zoom + value < 1:
            zoom = 1
        else:
            zoom += value
        resetProjectionMatrix(window.width)
    def apply_to_face_index(value):
        global flags, draw_mesh_index, draw_face_index, meshes
        if flags['mesh_drawn_since_last_face_index_change']:
            flags['draw_mesh'] = True
            number_of_faces = len(meshes[draw_mesh_index].face_meshes)
            if draw_face_index + value >= 0 and draw_face_index + value <= number_of_faces:
                draw_face_index += value
                flags['mesh_drawn_since_last_face_index_change'] = False
    if motion == key.MOTION_UP:
        apply_to_zoom(-ZOOM_STEP)
    if motion == key.MOTION_DOWN:
        apply_to_zoom(ZOOM_STEP)
    if motion == key.MOTION_LEFT:
        apply_to_rotation(-ROTATION_STEP_FAST)
    if motion == key.MOTION_PREVIOUS_WORD:
        apply_to_rotation(-ROTATION_STEP_SLOW)
    if motion == key.MOTION_RIGHT:
        apply_to_rotation(ROTATION_STEP_FAST)
    if motion == key.MOTION_NEXT_WORD:
        apply_to_rotation(ROTATION_STEP_SLOW)
    if motion == key.MOTION_PREVIOUS_PAGE:
        apply_to_face_index(1)
    if motion == key.MOTION_NEXT_PAGE:
        apply_to_face_index(-1)
    return pyglet.event.EVENT_HANDLED

@window.event
def on_key_press(symbol, modifiers):
    def reset():
        global zoom, rotation, tilt
        zoom = ZOOM_DEFAULT
        rotation = ROTATION_DEFAULT
        tilt = TILT_DEFAULT
    def apply_to_mesh_index(value):
        global draw_mesh_index, draw_face_index, meshes
        if draw_mesh_index + value >= 0 and draw_mesh_index + value < len(meshes):
            draw_mesh_index += value
            draw_face_index = 0
    global draw_face_index, rotation, tilt, flags
    if modifiers == 0:
        if symbol == key.C:
            reset()
        if symbol == key.A:
            flags['arrow_mode'] = not flags['arrow_mode']
        if symbol == key.K:
            apply_to_mesh_index(1)
        if symbol == key.J:
            apply_to_mesh_index(-1)
        if symbol == key.U:
            load_meshes_from_files()
        if symbol == key.V:
            flags['draw_vertices'] = not flags['draw_vertices']
        if symbol == key.B:
            flags['draw_bounding_box'] = not flags['draw_bounding_box']
        if symbol == key.N:
            flags['draw_mesh'] = not flags['draw_mesh']
            draw_face_index = 0
        if symbol == key.T:
            tilt = 10 if tilt == 0 else 0
    #return pyglet.event.EVENT_HANDLED # disables ESC termination handler

def update_direction_of_flight():
    global direction_of_flight
    matrix = (GLfloat * 16)()
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    x = matrix[2]
    y = matrix[6]
    z = matrix[10]
    direction_of_flight = np.array([x, y, z])

def normalize(vector, axis=-1, order=2):
    l2 = np.array(np.linalg.norm(vector, order, axis))
    l2[l2==0] = 1
    return vector / np.expand_dims(l2, axis)

@window.event
def on_draw():
    def set_modelview(mesh1D):
        global rotation, tilt
        scale = SAVE_ZONE_FACTOR * ZOOM_DEFAULT * 1.0/mesh1D.get_bb_size_factor()
        center = mesh1D.get_bb_center()
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(tilt, 1.0, 0.0, 0.0)
        glRotatef(rotation, 0.0, 1.0, 0.0)
        glScalef(scale, scale, scale)
        glTranslatef(-center[0], -center[1], -center[2])
        update_direction_of_flight()
    def draw_bb(mesh1D):
        bounding_box = mesh1D.bounding_box
        x_min, x_max, y_min, y_max, z_min, z_max = bounding_box
        glColor3f(0.5, 0.5, 0.5)
        glBegin(GL_LINE_LOOP)
        glVertex3f(x_min, y_min, z_min)
        glVertex3f(x_max, y_min, z_min)
        glVertex3f(x_max, y_max, z_min)
        glVertex3f(x_min, y_max, z_min)
        glEnd()
        glBegin(GL_LINE_LOOP)
        glVertex3f(x_min, y_min, z_max)
        glVertex3f(x_max, y_min, z_max)
        glVertex3f(x_max, y_max, z_max)
        glVertex3f(x_min, y_max, z_max)
        glEnd()
        glBegin(GL_LINES)
        glVertex3f(x_min, y_min, z_min)
        glVertex3f(x_min, y_min, z_max)
        glVertex3f(x_max, y_min, z_min)
        glVertex3f(x_max, y_min, z_max)
        glVertex3f(x_max, y_max, z_min)
        glVertex3f(x_max, y_max, z_max)
        glVertex3f(x_min, y_max, z_min)
        glVertex3f(x_min, y_max, z_max)
        glEnd()
    def draw_mesh(mesh1D):
        def draw_wire_mesh(vertices, wire_mesh):
            def draw_line(v0, v1):
                global flags
                if flags['arrow_mode'] and draw_face_index > 0:
                    global direction_of_flight
                    v0 = np.array(v0)
                    v1 = np.array(v1)
                    v01 = v1-v0
                    v_dof = normalize(direction_of_flight)
                    v_norm = np.cross(v01, v_dof)
                    width_factor = 0.05
                    length_factor = 0.4
                    va = v0 + (1-length_factor)*v01 + width_factor * v_norm
                    vb = v0 + (1-length_factor)*v01 + width_factor * -v_norm
                    # draw line
                    glVertex3f(v0[0], v0[1], v0[2])
                    glVertex3f(v1[0], v1[1], v1[2])
                    # draw arrow wing a
                    glVertex3f(va[0], va[1], va[2])
                    glVertex3f(v1[0], v1[1], v1[2])
                    # draw arrow wing b
                    glVertex3f(vb[0], vb[1], vb[2])
                    glVertex3f(v1[0], v1[1], v1[2])
                else:
                    glVertex3f(v0[0], v0[1], v0[2])
                    glVertex3f(v1[0], v1[1], v1[2])
            glBegin(GL_LINES)
            for iVertex in range(0, len(wire_mesh)):
                v0 = vertices[wire_mesh[iVertex]]
                v1 = vertices[wire_mesh[(iVertex+1)%len(wire_mesh)]]
                draw_line(v0, v1)
            glEnd()
        vertices = mesh1D.vertices
        face_meshes = mesh1D.face_meshes
        glColor3f(1.0, 1.0, 1.0)
        for iFace in range(0, len(face_meshes)):
            if draw_face_index > 0 and draw_face_index != iFace:
                continue
            face_mesh = face_meshes[iFace]
            outer_wire_mesh, inner_wire_meshes = face_mesh
            draw_wire_mesh(vertices, outer_wire_mesh)
            for inner_wire_mesh in inner_wire_meshes:
                draw_wire_mesh(vertices, inner_wire_mesh)
        flags['mesh_drawn_since_last_face_index_change'] = True
    def draw_vertices(mesh1D):
        glColor3f(1.0, 0.0, 0.0)
        glBegin(GL_POINTS)
        for v in mesh1D.vertices:
            x, y, z = v
            glVertex3f(x, y, z)
        glEnd()
    def draw_label(label_id, text):
        font_size = 14
        padding = 10

        x_min = -window.width//2
        x_pos = x_min + padding

        y_max = window.height//2
        y_pos = (y_max - padding) - label_id * (font_size + padding)

        label = pyglet.text.Label(text, font_name='Arial', font_size=font_size,
                                  x=x_pos, y=y_pos, anchor_x='left', anchor_y='top')
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        resetProjectionMatrix(window.width)
        label.draw()
    global flags, zoom, meshes, draw_mesh_index
    window.set_caption(window_caption)
    if len(meshes) == 0:
        return
    mesh1D = meshes[draw_mesh_index]
    glClear(GL_COLOR_BUFFER_BIT)
    resetProjectionMatrix(zoom)
    set_modelview(mesh1D)
    if flags['draw_bounding_box']:
        draw_bb(mesh1D)
    if flags['draw_mesh'] and draw_face_index == 0:
        draw_mesh(mesh1D)
    if flags['draw_vertices']:
        draw_vertices(mesh1D)
    if flags['draw_mesh'] and draw_face_index != 0:
        draw_mesh(mesh1D)
    draw_label(0, mesh1D.name)
    draw_label(1, 'face index: '+str(draw_face_index))

if __name__ == '__main__':
    load_meshes_from_files()
    pyglet.app.run()
