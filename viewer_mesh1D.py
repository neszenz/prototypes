import os
import glob
import numpy as np
import pickle
import sys

import pyglet
from pyglet.gl import *
from pyglet.window import key, mouse

import mesh1D
from mesh1D import SuperVertex, Mesh1D

## config constants  = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_DIR = 'tmp'
DEFAULT_X = 1280
DEFAULT_Y = 720
ZOOM_DEFAULT = 100
ZOOM_STEP = ZOOM_DEFAULT / 10
SAVE_ZONE_FACTOR = 0.8 # scale weight for model size to window ratio
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
gvars = {
    'mesh_index' : 0,
    'face_index' : 0, # 0: all faces; >1: select individual face
    'dof' : np.array([0.0, 0.0, 1.0]), # direction of flight
    'zoom' : ZOOM_DEFAULT,
    'rotation' : np.array([0.0, 0.0, 0.0]),
    'tilt' : TILT_DEFAULT,
    'drag' : np.array([0.0, 0.0, 0.0])
}

def load_meshes_from_files():
    global meshes
    meshes.clear()
    files = sorted(glob.glob(os.path.join(INPUT_DIR, '**/*'+mesh1D.FILE_EXTENSION)), reverse=True)
    for path in files:
        loaded_mesh1D = pickle.load(open(path, 'rb'))
        loaded_mesh1D.name = path
        meshes.append(loaded_mesh1D)

def set_projection(width):
    w = width
    h = w / (window.width/window.height)
    d = 99999.0
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(-w/2, w/2, -h/2, h/2, -d/2, d/2)
    glMatrixMode(GL_MODELVIEW)

@window.event
def on_resize(width, height):
    glViewport(0, 0, width, height)
    return pyglet.event.EVENT_HANDLED

def apply_to_rotation(value, axis):
    global gvars
    if gvars['rotation'][axis] + value < 0:
        gvars['rotation'][axis] = 360 + value
    elif gvars['rotation'][axis] + value > 360:
        gvars['rotation'][axis] = value
    elif gvars['rotation'][axis] + value == 360:
            gvars['rotation'][axis] = 0
    else:
        gvars['rotation'][axis] += value
def apply_to_zoom(value):
    global gvars
    if gvars['zoom'] + value < 1:
        gvars['zoom'] = 1
    else:
        gvars['zoom'] += value

@window.event
def on_text_motion(motion):
    def apply_to_face_index(value):
        global flags, gvars, meshes
        if flags['mesh_drawn_since_last_face_index_change']:
            flags['draw_mesh'] = True
            number_of_faces = len(meshes[gvars['mesh_index']].face_meshes)
            if gvars['face_index'] + value >= 0 and gvars['face_index'] + value <= number_of_faces:
                gvars['face_index'] += value
                flags['mesh_drawn_since_last_face_index_change'] = False
    if motion == key.MOTION_UP:
        apply_to_zoom(-ZOOM_STEP)
    if motion == key.MOTION_DOWN:
        apply_to_zoom(ZOOM_STEP)
    if motion == key.MOTION_LEFT:
        apply_to_rotation(-ROTATION_STEP_FAST, 1)
    if motion == key.MOTION_PREVIOUS_WORD:
        apply_to_rotation(-ROTATION_STEP_SLOW, 1)
    if motion == key.MOTION_RIGHT:
        apply_to_rotation(ROTATION_STEP_FAST, 1)
    if motion == key.MOTION_NEXT_WORD:
        apply_to_rotation(ROTATION_STEP_SLOW, 1)
    if motion == key.MOTION_PREVIOUS_PAGE:
        apply_to_face_index(1)
    if motion == key.MOTION_NEXT_PAGE:
        apply_to_face_index(-1)
    return pyglet.event.EVENT_HANDLED

@window.event
def on_key_press(symbol, modifiers):
    def reset():
        global gvars
        if gvars['zoom'] != ZOOM_DEFAULT or np.linalg.norm(gvars['drag']) > 0:
            gvars['zoom'] = ZOOM_DEFAULT
            gvars['drag'] = np.array([0.0, 0.0, 0.0])
        else:
            gvars['rotation'] = np.array([0.0, 0.0, 0.0])
            gvars['tilt'] = TILT_DEFAULT
    def apply_to_mesh_index(value):
        global gvars, meshes
        if gvars['mesh_index'] + value >= 0 and gvars['mesh_index']+ value < len(meshes):
            gvars['mesh_index']+= value
    global flags, gvars
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
            gvars['face_index'] = 0
        if symbol == key.T:
            gvars['tilt'] = 10 if gvars['tilt'] == 0 else 0
        if symbol == key.X:
            pass
    #return pyglet.event.EVENT_HANDLED # disables ESC termination handler

@window.event
def on_mouse_drag(x, y, dx, dy, button, modifiers):
    global gvars
    scale = gvars['zoom']/window.width
    if button == mouse.LEFT and modifiers == 0:
        if gvars['rotation'][0] > 90 and gvars['rotation'][0] < 270:
            apply_to_rotation(-1*dx*scale, 1)
        else:
            apply_to_rotation(dx*scale, 1)
        apply_to_rotation(-1*dy*scale, 0)
    if button == mouse.RIGHT or (button == mouse.LEFT and modifiers == key.MOD_ALT):
        gvars['drag'][0] += dx * scale
        gvars['drag'][1] += dy * scale

@window.event
def on_mouse_scroll(x, y, scroll_x, scroll_y):
    apply_to_zoom(-3*scroll_y)

def update_direction_of_flight():
    global gvars
    matrix = (GLfloat * 16)()
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    x = matrix[2]
    y = matrix[6]
    z = matrix[10]
    gvars['dof'] = np.array([x, y, z])

def normalize(vector, axis=-1, order=2):
    l2 = np.array(np.linalg.norm(vector, order, axis))
    l2[l2==0] = 1
    return vector / np.expand_dims(l2, axis)

@window.event
def on_draw():
    def set_modelview(mesh1D):
        global gvars
        scale = SAVE_ZONE_FACTOR * ZOOM_DEFAULT * 1.0/mesh1D.get_bb_size_factor()
        drag = gvars['drag']
        center = mesh1D.get_bb_center()
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glTranslatef(drag[0], drag[1], drag[2])
        glScalef(scale, scale, scale)
        glRotatef(gvars['tilt'], 1.0, 0.0, 0.0)
        glRotatef(gvars['rotation'][0], 1.0, 0.0, 0.0)
        glRotatef(gvars['rotation'][1], 0.0, 1.0, 0.0)
        glRotatef(gvars['rotation'][2], 0.0, 0.0, 1.0)
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
            def draw_line(sv0, sv1):
                global flags
                v0 = sv0.to_vec3()
                v1 = sv1.to_vec3()
                # draw line segment
                glVertex3f(v0[0], v0[1], v0[2])
                glVertex3f(v1[0], v1[1], v1[2])
                if flags['arrow_mode'] and gvars['face_index'] > 0:
                    v01 = v1-v0
                    v_dof = normalize(gvars['dof'])
                    v_norm = np.cross(v01, v_dof)
                    width_factor = 0.05
                    length_factor = 0.4
                    va = v0 + (1-length_factor)*v01 + width_factor * v_norm
                    vb = v0 + (1-length_factor)*v01 + width_factor * -v_norm
                    # draw arrow wing a and b
                    glVertex3f(va[0], va[1], va[2])
                    glVertex3f(v1[0], v1[1], v1[2])
                    glVertex3f(vb[0], vb[1], vb[2])
                    glVertex3f(v1[0], v1[1], v1[2])
            glBegin(GL_LINES)
            for iVertex in range(0, len(wire_mesh)):
                sv0 = vertices[wire_mesh[iVertex]]
                sv1 = vertices[wire_mesh[(iVertex+1)%len(wire_mesh)]]
                draw_line(sv0, sv1)
            glEnd()
        vertices = mesh1D.vertices
        face_meshes = mesh1D.face_meshes
        glColor3f(1.0, 1.0, 1.0)
        for iFace in range(0, len(face_meshes)):
            if gvars['face_index'] > 0 and gvars['face_index']-1 != iFace:
                continue
            face_mesh = face_meshes[iFace]
            outer_wire_mesh, inner_wire_meshes, _ = face_mesh
            draw_wire_mesh(vertices, outer_wire_mesh)
            for inner_wire_mesh in inner_wire_meshes:
                draw_wire_mesh(vertices, inner_wire_mesh)
        flags['mesh_drawn_since_last_face_index_change'] = True
    def draw_vertices(mesh1D):
        glColor3f(1.0, 0.0, 0.0)
        glBegin(GL_POINTS)
        for sv in mesh1D.vertices:
            x, y, z = sv.to_vec3()
            glVertex3f(x, y, z)
        glEnd()
    def draw_labels(mesh):
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
            set_projection(window.width)
            label.draw()
        global gvars
        draw_label(0, 'draw mesh '+str(-1*gvars['mesh_index'])+' ('+mesh.name+')')
        if gvars['face_index'] == 0:
            draw_label(1, 'draw all faces')
        else:
            face_id_from_total = str(gvars['face_index']) + '/' + str(len(mesh.face_meshes))
            face_type = mesh.get_face_type(gvars['face_index']-1)
            draw_label(1, 'draw face '+face_id_from_total+' ('+face_type+')')
    global flags, gvars
    window.set_caption(window_caption)
    if len(meshes) == 0:
        return
    mesh1D = meshes[gvars['mesh_index']]
    glClear(GL_COLOR_BUFFER_BIT)
    set_projection(gvars['zoom'])
    set_modelview(mesh1D)
    if flags['draw_bounding_box']:
        draw_bb(mesh1D)
    if flags['draw_mesh'] and gvars['face_index'] == 0:
        draw_mesh(mesh1D)
    if flags['draw_vertices']:
        draw_vertices(mesh1D)
    if flags['draw_mesh'] and gvars['face_index'] != 0:
        draw_mesh(mesh1D)
    draw_labels(mesh1D)

if __name__ == '__main__':
    load_meshes_from_files()
    pyglet.app.run()
