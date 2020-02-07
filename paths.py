import os
import pathlib

RESOURCE_DIR = '../resources'
RESULTS_DIR = '../results'

TMP_DIR = os.path.join(RESULTS_DIR, 'tmp')

# original data set collection
BOX = os.path.join(RESOURCE_DIR, 'step/boxWithHole.step')
CYLINDER = os.path.join(RESOURCE_DIR, 'step/cylinder.step')
FILLET = os.path.join(RESOURCE_DIR, 'step/fillet.step')
NUT = os.path.join(RESOURCE_DIR, 'step/hybrid_occ_builtin.step')
PIPES = os.path.join(RESOURCE_DIR, 'step/pipes.step')
REVOLVE = os.path.join(RESOURCE_DIR, 'step/revolve2.step')
SPHERE = os.path.join(RESOURCE_DIR, 'step/spherical_surf.step')
SURFFIL = os.path.join(RESOURCE_DIR, 'step/surface_filling.step')
TWIST = os.path.join(RESOURCE_DIR, 'step/twist.step')
TEST1 = os.path.join(RESOURCE_DIR, 'step/zTest1.step')

# corrupt data set
#STEP_X = '../resources/corrupt/00109971_343538c408e4f4299b3e20a4_step_001.step'
X = '../resources/corrupt/00100003_35809b42367c59ffb8aca25e_step_001.step'

# custom files created with gmsh
def custom(index):
    assert index >= 0

    base_path = os.path.join(RESOURCE_DIR, 'gmsh_custom/')
    file_name = 'c'+str(index)+'.step'

    paths = list(pathlib.Path(base_path).rglob(file_name))
    assert len(paths) == 1

    return str(paths[0])

# abc data set collection: 42, 111, 9999, 8613
def abc(index):
    assert index >= 0
    index_str = str(index)
    assert len(index_str) <= 4

    base_path = os.path.join(RESOURCE_DIR, 'abc/step/')
    padding = '000'[0:4-len(index_str)]
    regex_name = '0009' + padding + index_str + '/*.step'

    paths = list(pathlib.Path(base_path).rglob(regex_name))
    assert not len(paths) == 0
    assert not len(paths) > 1

    return str(paths[0])
