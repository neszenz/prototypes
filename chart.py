import matplotlib.pyplot as plt
import numpy as np

classes = (
    'Plane',
    'Cylinder',
    'Cone',
    'Shpere',
    'Torus',
    'Bézier',
    'B-Spline',
    'Revolution',
    'Extrusion',
    'Offset',
    'Other'
)
y_pos = np.arange(len(classes))
amounts = [
    246962,
    142907,
    12505,
    5325,
    21616,
    0,
    40387,
    2213,
    34808,
    966,
    0
]

plt.bar(y_pos, amounts, align='center', alpha=0.8)
plt.xticks(y_pos, classes, rotation='vertical')
plt.ylabel('number of occurrence')
plt.yscale('linear')
plt.title('surface type distribution of 10th abc data chunk (w/o #91715)')

plt.show()
