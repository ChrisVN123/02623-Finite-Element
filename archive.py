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


def u_fem_global_sum_exp(xq, nodes):
    """
    Interpolant of exp(x) on the mesh given by `nodes`:
      u_h(x) = sum_i exp(x_i) * N_i(x)
    """
    nodes = np.asarray(nodes, dtype=float)
    coeffs = np.exp(nodes)

    xq_arr = np.asarray(xq, dtype=float)
    uh = np.zeros_like(xq_arr, dtype=float)

    for i in range(len(nodes)):
        uh += coeffs[i] * N_hat(xq_arr, i, nodes)

    return float(uh) if np.isscalar(xq) else uh

@nb.njit(parallel=True, fastmath=True)
def create_A_diagonals(h_e):
    ne = h_e.size
    n = ne + 1

    diag  = np.empty(n, dtype=np.float64)
    upper = np.empty(ne, dtype=np.float64) 
    lower = np.empty(ne, dtype=np.float64)  

    for e in nb.prange(ne):
        h = h_e[e]
        v = (-1.0 / h) + (h / 6.0)
        upper[e] = v
        lower[e] = v

    diag[0] = (1.0 / h_e[0]) + (h_e[0] / 3.0)
    diag[n - 1] = (1.0 / h_e[ne - 1]) + (h_e[ne - 1] / 3.0)

    for k in nb.prange(1, n - 1):
        hl = h_e[k - 1]
        hr = h_e[k]
        diag[k] = (1.0 / hl + hl / 3.0) + (1.0 / hr + hr / 3.0)

    return lower, diag, upper

def create_A(h_e):
    h_e = np.asarray(h_e, dtype=np.float64)
    lower, diag, upper = create_A_diagonals(h_e)
    return diags([lower, diag, upper], offsets=[-1, 0, 1], format="csr")
