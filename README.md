# Prototypes
The various Python modules in this repository are partly dependent and partly independent on one another. They contain all the Python code I program during my thesis on Delaunay Refinement of smooth surface domains embedded into 3-dimensional space.

Dependencies:
- [pythonOCC](https://github.com/tpaviot/pythonocc)
- [pytriangle](https://github.com/pletzer/pytriangle)
- [pyglet](http://pyglet.org)
- [OpenMesh](https://openmesh-python.readthedocs.io/en/latest/)
- [gerobust](https://github.com/aluriak/gerobust)

## Results

Recreation of figure from Chew's 1993s paper for the 2D algorithm (dr2d.py).
![](imgs/chew93_2D.png "recreation of chew93 figure for planar domain")

Results of my implementation of surface Delaunay refinement (mesher_v2.py) based on Chew's paper.
![](imgs/chew93_Surface.png "various results from mesher_v2.py")

After the initial triangulation (A) already approximates geometry and topology satisfactorialy, iteratively refining it leads to fold generation (B) and in this specific case even to an infinite loop in the travel test (C).
![](imgs/fold_generation.png "one of the error cases lead to infinite looping")

The range of mesh quality (measured as radius-edge-ratio) convergest during refinement to the interval [1/sqrt(3), 1].
![](imgs/radius_edge_ratio.png "mesh quality convergence")

The size threshold helps to guarantee interior features to be sampled appropriately.
![](imgs/size_test.png "results for various size thresholds")

CAD surface type distribution (filter_v2.py, char.py) illustrates benefit of switching between 2D and surface algorithm.
![](imgs/surface_type_distribution_chart.png "CAD surface type distribution")
