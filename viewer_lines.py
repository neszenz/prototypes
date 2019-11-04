import os
import glob
import numpy as np

import pyglet
from pyglet.gl import *
from pyglet.window import key
from pyglet.window import mouse

INPUT_DIR = 'tmp'

DEFAULT_X = 1280
DEFAULT_Y = 720
zoom = 2
ZOOM_STEP = 0.1
rotation = 0
ROTATION_STEP = 10
draw_bounding_box = False

lines = list()
bounding_box = (
    float('inf'), float('-inf'),
    float('inf'), float('-inf'),
    float('inf'), float('-inf')
)
window = pyglet.window.Window(DEFAULT_X, DEFAULT_Y,
                              fullscreen=False,
                              resizable=True,
                              vsync=True)
window_caption = 'Delaunay Surface Mesh Refinement'


def get_bb_size():
    global bounding_box
    x_min, x_max, y_min, y_max, z_min, z_max = bounding_box
    return (x_max-x_min, y_max-y_min, z_max-z_min)

def get_bb_size_factor():
    return np.linalg.norm(get_bb_size())

def get_bb_center():
    x_min, x_max, y_min, y_max, z_min, z_max = bounding_box
    cog = np.array((x_max-x_min, y_max-y_min, z_max-z_min)) / 2
    return np.array((x_min, y_min, z_min)) + cog

def get_path_to_latest_lines_file():
    files = sorted(glob.glob(os.path.join(INPUT_DIR, '**/*.lines')))
    return files.pop()

# lines are a tuple of two points, which are lists of coordinate values
# returns a list of those lines
def load_lines_form_file():
    def clear_lines():
        global lines, bounding_box
        lines.clear()
        bounding_box = (
            float('inf'), float('-inf'),
            float('inf'), float('-inf'),
            float('inf'), float('-inf')
        )
    def update_bounding_box_size(p):
        global bounding_box
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
        bounding_box = x_min, x_max, y_min, y_max, z_min, z_max
    global lines
    clear_lines()
    path = get_path_to_latest_lines_file()
    with open(path, 'r') as f:
        for l in f:
            fl = [float(i) for i in l.split()]
            update_bounding_box_size(fl[0:3])
            update_bounding_box_size(fl[3:])
            line = (fl[0:3], fl[3:])
            lines.append(line)

def updateProjectionMatrix(width, height):
    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    global zoom
    w = zoom
    h = w / (width/height)
    glOrtho(-w/2, w/2, -h/2, h/2, -w, w)
    glMatrixMode(GL_MODELVIEW)

@window.event
def on_resize(width, height):
    updateProjectionMatrix(width, height)
    return pyglet.event.EVENT_HANDLED

@window.event
def on_key_press(symbol, modifiers):
    def apply_to_zoom(value):
        global zoom
        if zoom + value < 1:
            zoom = 1
        else:
            zoom += value
        updateProjectionMatrix(window.width, window.height)
    global zoom, rotation, draw_bounding_box
    if symbol == key.UP:
        apply_to_zoom(-ZOOM_STEP)
    if symbol == key.DOWN:
        apply_to_zoom(ZOOM_STEP)
    if symbol == key.LEFT:
        rotation -= ROTATION_STEP
    if symbol == key.RIGHT:
        rotation += ROTATION_STEP
    if symbol == key.U and modifiers == 0:
        load_lines_form_file()
    if symbol == key.B and modifiers == 0:
        draw_bounding_box = not draw_bounding_box
    return pyglet.event.EVENT_HANDLED # disables ESC termination handler

@window.event
def on_draw():
    def set_modelview():
        global rotation
        scale_factor = 1.0/get_bb_size_factor()
        center = get_bb_center()
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glColor3f(1.0, 1.0, 1.0)
        glRotatef(10.0, 1.0, 0.0, 0.0)
        glRotatef(rotation, 0.0, 1.0, 0.0)
        glScalef(scale_factor, scale_factor, scale_factor)
        glTranslatef(-center[0], -center[1], -center[2])
        pass
    def draw_lines():
        glBegin(GL_LINES)
        for line in lines:
            x0, y0, z0 = line[0]
            x1, y1, z1 = line[1]
            glVertex3f(x0, y0, z0)
            glVertex3f(x1, y1, z1)
        glEnd()
    def draw_bb():
        global bounding_box
        x_min, x_max, y_min, y_max, z_min, z_max = bounding_box
        glColor3f(0.5, 0.0, 0.0)
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
    global draw_bounding_box
    glClear(GL_COLOR_BUFFER_BIT)
    set_modelview()
    draw_lines()
    if draw_bounding_box:
        draw_bb()

if __name__ == '__main__':
    load_lines_form_file()
    pyglet.app.run()
