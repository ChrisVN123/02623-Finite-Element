import numpy as np

h = 3
L = 2
c = 1
d = 2.71**2
M = 3

K = np.array([[1/h+h/3, -1/h+h/6],[-1/h+h/6, 1/h+h/3]])

A = np.zeros((3,3))
print()