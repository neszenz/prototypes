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
ZOOM_DEFAULT = 1.3
ZOOM_STEP = 0.1
ROTATION_DEFAULT = 0
ROTATION_STEP = 10
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
    'draw_lines' : True
}
draw_mesh_index = 0
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

@window.event
def on_key_press(symbol, modifiers):
    def reset():
        global zoom, rotation, tilt
        zoom = ZOOM_DEFAULT
        rotation = ROTATION_DEFAULT
        tilt = TILT_DEFAULT
    def apply_to_zoom(value):
        global zoom
        if zoom + value < 1:
            zoom = 1
        else:
            zoom += value
        resetProjectionMatrix(window.width)
    global draw_mesh_index, zoom, rotation, tilt, flags
    if symbol == key.J and modifiers == 0:
        draw_mesh_index -= 1 if draw_mesh_index > 0 else 0
    if symbol == key.K and modifiers == 0:
        draw_mesh_index += 1 if draw_mesh_index+1 < len(meshes) else 0
    if symbol == key.C and modifiers == 0:
        reset()
    if symbol == key.UP:
        apply_to_zoom(-ZOOM_STEP)
    if symbol == key.DOWN:
        apply_to_zoom(ZOOM_STEP)
    if symbol == key.LEFT:
        rotation -= ROTATION_STEP
    if symbol == key.RIGHT:
        rotation += ROTATION_STEP
    if symbol == key.U and modifiers == 0:
        load_meshes_from_files()
    if symbol == key.V and modifiers == 0:
        flags['draw_vertices'] = not flags['draw_vertices']
    if symbol == key.B and modifiers == 0:
        flags['draw_bounding_box'] = not flags['draw_bounding_box']
    if symbol == key.N and modifiers == 0:
        flags['draw_lines'] = not flags['draw_lines']
    if symbol == key.T and modifiers == 0:
        tilt = 10 if tilt == 0 else 0
    return pyglet.event.EVENT_HANDLED # disables ESC termination handler

@window.event
def on_draw():
    def set_modelview(mesh1D):
        global rotation, tilt
        scale_factor = 1.0/mesh1D.get_bb_size_factor()
        center = mesh1D.get_bb_center()
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(tilt, 1.0, 0.0, 0.0)
        glRotatef(rotation, 0.0, 1.0, 0.0)
        glScalef(scale_factor, scale_factor, scale_factor)
        glTranslatef(-center[0], -center[1], -center[2])
        pass
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
            glBegin(GL_LINE_LOOP)
            for vi in wire_mesh:
                x, y, z = vertices[vi]
                glVertex3f(x, y, z)
            glEnd()
        vertices = mesh1D.vertices
        face_meshes = mesh1D.face_meshes
        glColor3f(1.0, 1.0, 1.0)
        for face_mesh in face_meshes:
            outer_wire_mesh, inner_wire_meshes = face_mesh
            draw_wire_mesh(vertices, outer_wire_mesh)
            for inner_wire_mesh in inner_wire_meshes:
                draw_wire_mesh(vertices, inner_wire_mesh)
    def draw_vertices(mesh1D):
        glColor3f(1.0, 0.0, 0.0)
        glBegin(GL_POINTS)
        for v in mesh1D.vertices:
            x, y, z = v
            glVertex3f(x, y, z)
        glEnd()
    def draw_label(text):
        label = pyglet.text.Label(text,
                                  font_name='Arial', font_size=14,
                                  x=10-window.width//2, y=-10+window.height//2,
                                  anchor_x='left', anchor_y='top')
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
    if flags['draw_lines']:
        draw_mesh(mesh1D)
    if flags['draw_vertices']:
        draw_vertices(mesh1D)
    draw_label(mesh1D.name)

if __name__ == '__main__':
    load_meshes_from_files()
    pyglet.app.run()
