import numpy as np
import numba as nb
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, diags
from scipy.sparse.linalg import spsolve
import time

L = 2.0
c = 1.0
d_bc = np.exp(2)



def create_A(h_e):
    ne = len(h_e)
    n_nodes = ne + 1
    A = lil_matrix((n_nodes, n_nodes), dtype=float)

    for e in range(ne):
        h = h_e[e]
        Ae = np.array([[ 1/h + h/3,  -1/h + h/6],
                       [ -1/h + h/6,  1/h + h/3]], dtype=float)

        i = e
        j = e + 1

        A[i, i] += Ae[0, 0]
        A[i, j] += Ae[0, 1]
        A[j, i] += Ae[1, 0]
        A[j, j] += Ae[1, 1]

    return A.tocsr()


def create_b(n_nodes):
    return np.zeros(n_nodes, dtype=float)

def gaussian_solve(A, b, left_bc, right_bc):
    n_nodes = A.shape[0]
    interior = np.arange(1, n_nodes - 1)

    # we now remove the columns and rows related to u_1 and u_N
    A_ii = A[interior][:, interior]
    b_i = b[interior].copy()

    # Dirichlet elimination
    b_i -= A[interior, 0].toarray().ravel()  * left_bc
    b_i -= A[interior, -1].toarray().ravel() * right_bc
    u_i = spsolve(A_ii, b_i)

    u = np.zeros(n_nodes, dtype=float)
    u[0] = left_bc
    u[-1] = right_bc
    u[interior] = u_i
    return u

def solve(ne, h=None):
    t0 = time.time()
    # create list of uniform interval lengths
    if h == None: h_e = np.full(ne, L / ne, dtype=float)
    else: 
        h_e = h
    # create list of x_values for plotting which is the accumulated sum of h_e
    x = np.r_[0.0, np.cumsum(h_e)]
    x[-1] = L  # ensure exact endpoint

    A = create_A(h_e)
    b = create_b(ne + 1)

    u = gaussian_solve(A, b, c, d_bc)
    t1 = time.time()
    return u, x, t1-t0

# ---- run ----

u_pde, x_nodes, time_g = solve(5)
print("PDE FEM nodal u:", u_pde)
print("nodes:", x_nodes)
print("Time spent:", time_g)

exit()
# testing runtime across iterations
n = 2000
runtime = np.zeros(n)
for i in range(n):
    _, _, time_g = solve(i+3)
    runtime[i]=time_g
xf = np.linspace(0, n, n)
plt.plot(xf,runtime)
plt.show()


# exact exp(x) for comparison with the interpolant
xf = np.linspace(0, L, 400)
exp_exact = np.exp(xf)

# interpolated exp(x) using basis on the FEM mesh x_nodes
uh_vals = u_fem_global_sum_exp(xf, x_nodes)

# plot
plt.figure()
plt.plot(x_nodes, u_pde, marker="o", label="FEM solution of PDE (nodal)")
plt.plot(xf, exp_exact, label="exp(x) exact")
plt.plot(xf, uh_vals, "--", label="exp(x) interpolated on mesh")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.grid(True)
plt.legend()
plt.show()

exit()
alpha = (d_bc - c*np.exp(-L)) / (np.exp(L) - np.exp(-L))
beta  = c - alpha

ne_values = np.arange(2, 1000)
h_values  = L / ne_values    

err_L2 = np.zeros_like(ne_values, dtype=float)
err_Linf = np.zeros_like(ne_values, dtype=float)

for k, ne in enumerate(ne_values):
    u_num, x_nodes = solve(ne)  

    u_ex = np.exp(x_nodes)
    e = u_num - u_ex

    err_L2[k] = np.sqrt(np.mean(e**2))

p_L2 = np.polyfit(np.log(h_values), np.log(err_L2), 1)[0]

print(f"Estimated convergence rate (RMS/L2): {p_L2:.3f}")

plt.figure()
plt.loglog(h_values, err_L2, marker="o", label="RMS (discrete L2) error")
plt.xlabel("h = L/ne")
plt.ylabel("error")
plt.grid(True, which="both")
plt.legend()
plt.show()
