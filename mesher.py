""" Chew93_Surface mesher for STEP files using sampler.py w/ OpenCASCADE
Contains functionality to generate a 2-dimensional mesh from a given path to a
STEP file.

To gain a foundation for the meshing process, sampler.py is used to generate a
1-dimensional mesh from the STEP file which is than iterated upon. Here, the
first step is to compute constrained Delaunay Triangulations w/ pytriangle and
then apply my implementation of Chew's Surface Delaunay Refinement algorithm.
"""
import triangle

import paths
import sampler

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_PATH = paths.PATH_98613
sampler.NUMBER_OF_SAMPLES = 20

if __name__ == '__main__':
    simple_sampler = sampler.factory(sampler.SAMPLER_TYPE.SIMPLE)
    mesh1D = simple_sampler(INPUT_PATH)
    sampler.write_mesh_to_file(mesh1D)
