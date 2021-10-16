import bezier
import numpy as np
from numpy import diff
import matplotlib.pyplot as plt


nodes = np.asfortranarray([
    [-1, -0.5, 0, 0.5, 1],
    [0, -0.5, 0, 0.5, 0],
])
curve = bezier.Curve(nodes, degree=4)

d = np.array([curve.evaluate_hodograph(t)[0:2] for t in np.arange(0, 1, 0.01)])

x = d[:, 0]
y = d[:, 1]

plt.plot(x, y)
plt.show()

