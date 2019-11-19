import os
import glob
import numpy as np
import pickle
import sys

import pyglet
from pyglet.gl import *
from pyglet.window import key, mouse

from meshkD import SuperVertex, MeshkD

## config constants  = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_DIR = 'tmp'
DEFAULT_X = 1280
DEFAULT_Y = 720
ZOOM_DEFAULT = 100
ZOOM_STEP = ZOOM_DEFAULT / 10
SAVE_ZONE_FACTOR = 0.7 # scale weight for model size to window ratio
ROTATION_STEP_SLOW = 1
ROTATION_STEP_FAST = 5
TILT_DEFAULT = 0

## global variables  = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
meshes = [] # list of tuples (mesh1D, mesh2D)
window = pyglet.window.Window(DEFAULT_X, DEFAULT_Y,
                              fullscreen=False,
                              resizable=True,
                              vsync=True)
window_caption = 'Delaunay Surface Mesh Refinement'

flags = {
    'draw_bounding_box' : True,
    'draw_vertices' : True,
    'draw_mesh1D' : True,
    'draw_mesh2D' : True,
    'draw_2d_mode' : False,
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

def num_of_faces_from_current_mesh():
    if len(meshes) == 0:
        return 0
    return meshes[gvars['mesh_index']].number_of_faces()

def load_meshes_from_files():
    global meshes
    meshes.clear()

    files = sorted(glob.glob(os.path.join(INPUT_DIR, '*'+MeshkD.FILE_EXTENSION)), reverse=True)

    for path in files:
        mesh = pickle.load(open(path, 'rb'))
        mesh.name = path + ' from ' + mesh.name

        meshes.append(mesh)

    gvars['face_index'] = min(gvars['face_index'], num_of_faces_from_current_mesh())

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
        global flags, gvars
        if flags['mesh_drawn_since_last_face_index_change']:
            number_of_faces = num_of_faces_from_current_mesh()
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
        if gvars['mesh_index']+value >= 0 and gvars['mesh_index']+value < len(meshes):
            gvars['mesh_index']+= value
        gvars['face_index'] = min(gvars['face_index'], num_of_faces_from_current_mesh())
    global flags, gvars
    if modifiers == 0:
        if symbol == key.C:
            reset()
        if symbol == key.A:
            flags['arrow_mode'] = not flags['arrow_mode']
        if symbol == key.F:
            flags['draw_2d_mode'] = not flags['draw_2d_mode']
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
            flags['draw_mesh1D'] = not flags['draw_mesh1D']
        if symbol == key.M:
            flags['draw_mesh2D'] = not flags['draw_mesh2D']
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
    def draw_no_meshes_msg():
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        set_projection(window.width)
        text = '<no meshes found in directory \'' + INPUT_DIR + '\'>'
        font_size = 14
        label = pyglet.text.Label(text, font_name='Arial', font_size=font_size,
                                  x=0.0, y=0.0, anchor_x='center', anchor_y='center')
        label.draw()
    def set_modelview(mesh):
        global gvars
        if flags['draw_2d_mode']:
            size_factor = mesh.get_bb2D_size_factor()
            center = mesh.get_bb2D_center()
        else:
            size_factor = mesh.get_bb3D_size_factor()
            center = mesh.get_bb3D_center()
        scale = SAVE_ZONE_FACTOR * ZOOM_DEFAULT * 1.0/size_factor
        drag = gvars['drag']
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
    def draw_bb(mesh):
        if flags['draw_2d_mode']:
            bounding_box = mesh.bounding_box2D
        else:
            bounding_box = mesh.bounding_box3D
        x_min, x_max, y_min, y_max, z_min, z_max = bounding_box
        glColor3f(0.5, 0.5, 0.0)
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
        def draw_face1D(face_mesh):
            def draw_wire(vertices, wire_mesh):
                def draw_line(sv0, sv1):
                    global flags
                    if flags['draw_2d_mode']:
                        v0 = sv0.UV_vec3()
                        v1 = sv1.UV_vec3()
                    else:
                        v0 = sv0.XYZ_vec3()
                        v1 = sv1.XYZ_vec3()
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
                for vertex_index in range(len(wire_mesh)):
                    sv0 = vertices[wire_mesh[vertex_index]]
                    sv1 = vertices[wire_mesh[(vertex_index+1)%len(wire_mesh)]]
                    draw_line(sv0, sv1)
                glEnd()
            vertices, wire_meshes, _, _ = face_mesh
            for wire_mesh in wire_meshes:
                draw_wire(vertices, wire_mesh)
        def draw_face2D(face_mesh):
            def draw_triangle(sv0, sv1, sv2):
                if flags['draw_2d_mode']:
                    v0 = sv0.UV_vec3()
                    v1 = sv1.UV_vec3()
                    v2 = sv2.UV_vec3()
                else:
                    v0 = sv0.XYZ_vec3()
                    v1 = sv1.XYZ_vec3()
                    v2 = sv2.XYZ_vec3()

                glVertex3f(v0[0], v0[1], v0[2])
                glVertex3f(v1[0], v1[1], v1[2])
                glVertex3f(v2[0], v2[1], v2[2])

                return
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
            vertices, _, triangles, _ = face_mesh

            glBegin(GL_TRIANGLES)

            for triangle in triangles:
                i0, i1, i2 = triangle
                draw_triangle(vertices[i0], vertices[i1], vertices[i2])

            glEnd()
            return
        for face_index in range(mesh.number_of_faces()):
            if gvars['face_index'] > 0 and gvars['face_index']-1 != face_index:
                continue
            face_mesh = mesh.face_meshes[face_index]
            if flags['draw_mesh2D']:
                glColor3f(0.6, 0.6, 0.6)
                draw_face2D(face_mesh)
            if flags['draw_mesh1D']:
                glColor3f(1.0, 1.0, 1.0)
                draw_face1D(face_mesh)
    def draw_vertices(mesh):
        def draw_super_vertex(sv):
            if flags['draw_2d_mode']:
                x, y, z = sv.UV_vec3()
            else:
                x, y, z = sv.XYZ_vec3()
            glVertex3f(x, y, z)

            return
        current_face_collection = []

        glBegin(GL_POINTS)

        for face_mesh in mesh.face_meshes:
            vertices, _, _, _ = face_mesh

            if gvars['face_index'] == 0:
                glColor3f(1.0, 1.0, 1.0)
            else:
                glColor3f(1.0, 0.0, 0.0)

            for sv in vertices:
                if gvars['face_index'] == sv.face_id:
                    current_face_collection.append(sv)
                draw_super_vertex(sv)

        glColor3f(1.0, 1.0, 1.0)
        for sv in current_face_collection:
            draw_super_vertex(sv)

        glEnd()
        return
    def draw_label_interface(mesh):
        def draw_label(label_id, text):
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            font_size = 14
            padding = 10

            x_min = -window.width//2
            x_pos = x_min + padding

            y_max = window.height//2
            y_pos = (y_max - padding) - label_id * (font_size + padding)

            label = pyglet.text.Label(text, font_name='Arial', font_size=font_size,
                                      x=x_pos, y=y_pos, anchor_x='left', anchor_y='top')
            label.draw()
        global gvars
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        set_projection(window.width)
        draw_label(0, 'mesh '+str(-1*gvars['mesh_index'])+' ('+mesh.name+')')
        if gvars['face_index'] == 0:
            draw_label(1, 'draw all faces')
        else:
            face_id_from_total = str(gvars['face_index']) + '/' + str(len(mesh.face_meshes))
            face_type = mesh.get_face_type(gvars['face_index']-1)
            draw_label(1, 'draw face '+face_id_from_total+' ('+face_type+')')
    global flags, gvars
    window.set_caption(window_caption)
    glClear(GL_COLOR_BUFFER_BIT)
    if len(meshes) == 0:
        draw_no_meshes_msg()
        return
    mesh = meshes[gvars['mesh_index']]
    set_projection(gvars['zoom'])
    set_modelview(mesh)
    if flags['draw_bounding_box']:
        draw_bb(mesh)
    if gvars['face_index'] == 0:
        draw_mesh(mesh)
    if flags['draw_vertices']:
        draw_vertices(mesh)
    if gvars['face_index'] != 0:
        draw_mesh(mesh)
    draw_label_interface(mesh)
    flags['mesh_drawn_since_last_face_index_change'] = True

if __name__ == '__main__':
    load_meshes_from_files()
    pyglet.app.run()
