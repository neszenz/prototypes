import os
import glob

import pyglet
from pyglet.gl import *
from pyglet.window import key
from pyglet.window import mouse

INPUT_DIR = 'tmp'

DEFAULT_X = 1280
DEFAULT_Y = 720
zoom = 25
ZOOM_STEP = 5
rotation = 0
ROTATION_STEP = 10

lines = list()
window = pyglet.window.Window(DEFAULT_X, DEFAULT_Y,
                              fullscreen=False,
                              resizable=True,
                              vsync=True)
window_caption = 'Delaunay Surface Mesh Refinement'


def get_path_to_latest_lines_file():
    files = sorted(glob.glob(os.path.join(INPUT_DIR, '**/*.lines')))
    return files.pop()

# lines are a tuple of two points, which are lists of coordinate values
# returns a list of those lines
def load_lines_form_file():
    global lines
    lines.clear()
    path = get_path_to_latest_lines_file()
    with open(path, 'r') as f:
        for l in f:
            fl = [float(i) for i in l.split()]
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
    global zoom, rotation
    if symbol == key.UP:
        apply_to_zoom(-ZOOM_STEP)
    if symbol == key.DOWN:
        apply_to_zoom(ZOOM_STEP)
    if symbol == key.LEFT:
        rotation -= ROTATION_STEP
    if symbol == key.RIGHT:
        rotation += ROTATION_STEP
    if symbol == key.U and modifiers == 0:
        global lines
        load_lines_form_file()
    # return pyglet.event.EVENT_HANDLED # disables ESC termination handler

@window.event
def on_draw():
    glClear(GL_COLOR_BUFFER_BIT)
    glLoadIdentity()
    global rotation
    glRotatef(10.0, 1.0, 0.0, 0.0)
    glRotatef(rotation, 0.0, 1.0, 0.0)
    glBegin(GL_LINES)
    for line in lines:
        x0, y0, z0 = line[0]
        x1, y1, z1 = line[1]
        glVertex3f(x0, y0, z0)
        glVertex3f(x1, y1, z1)
    glEnd()

if __name__ == '__main__':
    load_lines_form_file()
    pyglet.app.run()
