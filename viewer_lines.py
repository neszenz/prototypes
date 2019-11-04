import os
import glob
import numpy as np

import pyglet
from pyglet.gl import *
from pyglet.window import key
from pyglet.window import mouse

## config constants  = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_DIR = 'tmp'
DEFAULT_X = 1280
DEFAULT_Y = 720
ZOOM_DEFAULT = 1.1
ZOOM_STEP = 0.1
ROTATION_DEFAULT = 0
ROTATION_STEP = 10
TILT_DEFAULT = 0

## global variables  = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
meshes = list()
lines = list()
bounding_box = tuple()
window = pyglet.window.Window(DEFAULT_X, DEFAULT_Y,
                              fullscreen=False,
                              resizable=True,
                              vsync=True)
window_caption = 'Delaunay Surface Mesh Refinement'

flags = {
    'draw_bounding_box' : True,
    'draw_lines' : True
}
draw_mesh_index = 0
zoom = ZOOM_DEFAULT
rotation = ROTATION_DEFAULT
tilt = TILT_DEFAULT

class Mesh1D:
    def __init__(self, name):
        self.name = name
        self.lines = list()
        self.bounding_box = (
            float('inf'), float('-inf'),
            float('inf'), float('-inf'),
            float('inf'), float('-inf')
        )
    def get_bb_size(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box
        return (x_max-x_min, y_max-y_min, z_max-z_min)
    def get_bb_size_factor(self):
        return np.linalg.norm(self.get_bb_size())
    def get_bb_center(self):
        x_min, x_max, y_min, y_max, z_min, z_max = self.bounding_box
        cog = np.array((x_max-x_min, y_max-y_min, z_max-z_min)) / 2
        return np.array((x_min, y_min, z_min)) + cog

# lines are a tuple of two points, which are lists of coordinate values
# returns a list of those lines
def load_meshes_from_files():
    def load_mesh_from_file(path):
        def update_bounding_box(bounding_box, p):
            assert len(p) == 3
            x_min, x_max, y_min, y_max, z_min, z_max = bounding_box
            if p[0] < x_min:
                x_min = p[0]
            if p[0] > x_max:
                x_max = p[0]
            if p[1] < y_min:
                y_min = p[1]
            if p[1] > y_max:
                y_max = p[1]
            if p[2] < z_min:
                z_min = p[2]
            if p[2] > z_max:
                z_max = p[2]
            return x_min, x_max, y_min, y_max, z_min, z_max
        lines = list()
        bounding_box = (
            float('inf'), float('-inf'),
            float('inf'), float('-inf'),
            float('inf'), float('-inf')
        )
        with open(path, 'r') as f:
            for l in f:
                fl = [float(i) for i in l.split()]
                bounding_box = update_bounding_box(bounding_box, fl[0:3])
                bounding_box = update_bounding_box(bounding_box, fl[3:])
                line = (fl[0:3], fl[3:])
                lines.append(line)
        return lines, bounding_box
    global meshes
    meshes.clear()
    files = sorted(glob.glob(os.path.join(INPUT_DIR, '**/*.lines')), reverse=True)
    for path in files:
        new_mesh = Mesh1D(name=path)
        new_mesh.lines, new_mesh.bounding_box = load_mesh_from_file(path)
        meshes.append(new_mesh)

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
    if symbol == key.B and modifiers == 0:
        flags['draw_bounding_box'] = not flags['draw_bounding_box']
    if symbol == key.N and modifiers == 0:
        flags['draw_lines'] = not flags['draw_lines']
    if symbol == key.T and modifiers == 0:
        tilt = 10 if tilt == 0 else 0
    return pyglet.event.EVENT_HANDLED # disables ESC termination handler

@window.event
def on_draw():
    def set_modelview(mesh):
        global rotation, tilt
        scale_factor = 1.0/mesh.get_bb_size_factor()
        center = mesh.get_bb_center()
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(tilt, 1.0, 0.0, 0.0)
        glRotatef(rotation, 0.0, 1.0, 0.0)
        glScalef(scale_factor, scale_factor, scale_factor)
        glTranslatef(-center[0], -center[1], -center[2])
        pass
    def draw_bb(mesh):
        bounding_box = mesh.bounding_box
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
    def draw_mesh(mesh):
        lines = mesh.lines
        glColor3f(1.0, 1.0, 1.0)
        glBegin(GL_LINES)
        for line in lines:
            x0, y0, z0 = line[0]
            x1, y1, z1 = line[1]
            glVertex3f(x0, y0, z0)
            glVertex3f(x1, y1, z1)
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
    mesh = meshes[draw_mesh_index]
    glClear(GL_COLOR_BUFFER_BIT)
    resetProjectionMatrix(zoom)
    set_modelview(mesh)
    if flags['draw_bounding_box']:
        draw_bb(mesh)
    if flags['draw_lines']:
        draw_mesh(mesh)
    draw_label(mesh.name)

if __name__ == '__main__':
    load_meshes_from_files()
    pyglet.app.run()
