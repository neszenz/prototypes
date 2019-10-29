import pyglet
from pyglet.gl import *
from pyglet.window import key
from pyglet.window import mouse

DEFAULT_X = 1280
DEFAULT_Y = 720
PERSPECTIVE_MODE = False

window = pyglet.window.Window(DEFAULT_X, DEFAULT_Y,
                              fullscreen=False,
                              resizable=True,
                              vsync=True)
window_caption = 'Delaunay Surface Mesh Refinement'
# window.push_handlers(pyglet.window.event.WindowEventLogger()) # display events

@window.event
def on_draw():
    print('on_draw()')
    window.set_caption(window_caption)
    glClear(GL_COLOR_BUFFER_BIT)
    glLoadIdentity()
    vertices = (
        10.0, 10.0, 0.0,
        window.width-10.0, 10.0, 0.0,
        window.width-10.0, window.height*0.5, 0.0,
        window.width*0.5-5.0, window.height-10.0, 0.0,
        10.0, window.height*0.5, 0.0
    )
    colors = (
        1.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        1.0, 1.0, 1.0,
        0.0, 0.0, 1.0
    )
    normals = (
        0.0, 0.0, 1.0,
        0.0, 0.0, 1.0,
        0.0, 0.0, 1.0,
        0.0, 0.0, 1.0,
        0.0, 0.0, 1.0
    )
    indices = [0, 1, 2, 3, 4, 2, 0, 4, 1]
    dimensions = 3
    numOfVertices = int(len(vertices) / dimensions)
    vertex_data = ('v3f/dynamic', vertices)
    vertex_colors = ('c3f/static', colors)
    vertex_normals = ('n3f/stream', normals)
    vertex_list = pyglet.graphics.vertex_list_indexed(numOfVertices,
                                                      indices,
                                                      vertex_data,
                                                      vertex_colors,
                                                      vertex_normals)
    vertex_list.colors[9:12] = [1.0, 0.0, 1.0] # changing 4th vertex color
    vertex_list.draw(pyglet.gl.GL_LINE_LOOP)
    return pyglet.event.EVENT_HANDLED

def updateProjectionMatrix(width, height):
    if PERSPECTIVE_MODE:
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(65, width / float(height), .1, 1000)
        glMatrixMode(GL_MODELVIEW)
    else:
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(0, width, 0, height, -1, 1)
        glMatrixMode(GL_MODELVIEW)

@window.event
def on_resize(width, height):
    updateProjectionMatrix(width, height)
    return pyglet.event.EVENT_HANDLED

@window.event
def on_key_press(symbol, modifiers):
    print('on_key_press()')
    global window_caption
    if symbol == key.F and modifiers == 0:
        if window.fullscreen:
            window.set_fullscreen(False)
        else:
            window.set_fullscreen(True)
    if symbol == key.F and modifiers == 1:
        global PERSPECTIVE_MODE
        print(PERSPECTIVE_MODE)
        PERSPECTIVE_MODE = False if PERSPECTIVE_MODE else True
    # return pyglet.event.EVENT_HANDLED # disables ESC termination handler

@window.event
def on_mouse_press(x, y, button, modifiers):
    print('on_mouse_press()')
    return pyglet.event.EVENT_HANDLED

if __name__ == '__main__':
    pyglet.app.run()
    print('terminated!')
