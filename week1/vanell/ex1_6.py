import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return np.exp(-800*(x-0.4)**2)+0.25*np.exp(-40*(x-0.8)**2)


def N_hat(x, i, nodes):
    nodes = np.asarray(nodes, dtype=float)
    x_arr = np.asarray(x, dtype=float)
    M = len(nodes)

    if i < 0 or i >= M:
        raise ValueError("i out of range")
    if M < 2:
        raise ValueError("Need at least 2 nodes")

    if i == 0:
        x0, x1 = nodes[0], nodes[1]
        h = x1 - x0
        val = np.where((x_arr >= x0) & (x_arr <= x1), (x1 - x_arr) / h, 0.0)

    elif i == M - 1:
        x0, x1 = nodes[M - 2], nodes[M - 1]
        h = x1 - x0
        val = np.where((x_arr >= x0) & (x_arr <= x1), (x_arr - x0) / h, 0.0)

    else:
        xL, xC, xR = nodes[i - 1], nodes[i], nodes[i + 1]
        hL = xC - xL
        hR = xR - xC

        left  = np.where((x_arr >= xL) & (x_arr <= xC), (x_arr - xL) / hL, 0.0)
        right = np.where((x_arr >= xC) & (x_arr <= xR), (xR - x_arr) / hR, 0.0)
        val = left + right

    return float(val) if np.isscalar(x) else val

def x_interpolated(xq, nodes):
    """
    Interpolant of function(x) on the mesh given by `nodes`:
      u_h(x) = sum_i f(x_i) * N_i(x)
    """
    nodes = np.asarray(nodes, dtype=float)
    coeffs = func(nodes)

    xq_arr = np.asarray(xq, dtype=float)
    uh = np.zeros_like(xq_arr, dtype=float)

    for i in range(len(nodes)):
        uh += coeffs[i] * N_hat(xq_arr, i, nodes)

    return float(uh) if np.isscalar(xq) else uh

nodes = np.linspace(0,1,10)
x_val = np.linspace(0,1,100)


fx = x_interpolated(x_val,nodes)
print(fx)
func_vals = func(x_val)

plt.plot(x_val,func_vals)
plt.plot(x_val,fx)
plt.show()