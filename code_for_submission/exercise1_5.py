import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, diags
from scipy.sparse.linalg import spsolve
import time

L = 1.0
p = 1
epsilon = 0.01
start_bc = 0
end_bc = 0
f = 1

def create_A(h_e,eps):
    ne = len(h_e)
    n_nodes = ne + 1
    A = lil_matrix((n_nodes, n_nodes), dtype=float)


    for e in range(ne):
        h = h_e[e]
        Ae = np.array([[ p/2+eps/h, -eps/h+p/2 ],
                       [ -p/2-eps/h, eps/h-p/2]], dtype=float)

        i = e
        j = e + 1

        A[i, i] += Ae[0, 0]
        A[i, j] += Ae[0, 1]
        A[j, i] += Ae[1, 0]
        A[j, j] += Ae[1, 1]

    return A.tocsr()


def create_b(n_nodes):
    return np.zeros(n_nodes, dtype=float)

def gaussian_solve(A, b, left_bc, right_bc, h):
    n_nodes = A.shape[0]
    interior = np.arange(1, n_nodes - 1)

    # we now remove the columns and rows related to u_1 and u_N
    A_ii = A[interior][:, interior]
    b_i = b[interior].copy() + h

    b_i -= A[interior, 0].toarray().ravel()  * left_bc
    b_i -= A[interior, -1].toarray().ravel() * right_bc
    u_i = spsolve(A_ii, b_i)

    u = np.zeros(n_nodes, dtype=float)
    u[0] = left_bc
    u[-1] = right_bc
    u[interior] = u_i
    return u

def solve(ne, xn, bc0, bcn, eps, h=None):
    t0 = time.time()
    if h == None: h_e = np.full(ne, xn / ne, dtype=float)
    else: 
        h_e = h
    x = np.r_[0.0, np.cumsum(h_e)]
    x[-1] = xn  

    A = create_A(h_e,eps)
    b = create_b(ne + 1)

    u = gaussian_solve(A, b, bc0, bcn, h_e[0])
    t1 = time.time()
    return u, x, t1-t0

# ---- run ----

u_pde, x_nodes, time_g = solve(5000, xn=L, bc0=start_bc, bcn=end_bc, eps=epsilon)
print("PDE FEM nodal u:", u_pde)
print("nodes:", x_nodes)
print("Time spent:", time_g)


### d) 

u4, x4, _ = solve(1000, xn=L, bc0=start_bc, bcn=end_bc, eps=0.001)
u3, x3, _ = solve(100,  xn=L, bc0=start_bc, bcn=end_bc, eps=0.001)
u2, x2, _ = solve(50,   xn=L, bc0=start_bc, bcn=end_bc, eps=0.001)

plt.plot(x4, u4, label='ne=1000')
plt.plot(x3, u3, label='ne=100')
plt.plot(x2, u2, label='ne=50')
plt.legend()
plt.title('Approximating the solution using variable number of nodes')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.xlim(0, 1)
plt.grid(True, alpha=0.3)
plt.show()


# ## e) 

def function(x, phi, eps):
    return 1/phi*((1+(np.exp(phi/eps)-1))*x-np.exp(x*phi/eps))/(np.exp(phi/eps)-1)

n = 100
eps_list = [0.01]
error_con = np.zeros(n)
for e in eps_list:
    for i in range(3,n):
        u_fem,x_fem,_ = solve(i,xn=L, bc0=start_bc, bcn=end_bc, eps=e)
        u_func = function(x_fem,1,e)
        error_con[i] = np.max(np.abs(u_fem-u_func))
    plt.plot(error_con, label=f"eps={e:g}")
#plt.plot(error_con)
plt.title('Error analytical solution and FEM approximation')
plt.xlabel('Number of nodes')
plt.ylabel('Error')
plt.xlim(3,n)
plt.show()



# plot
u4, x4, _ = solve(20000, xn=L, bc0=start_bc, bcn=end_bc, eps=1)
u3, x3, _ = solve(20000,  xn=L, bc0=start_bc, bcn=end_bc, eps=0.01)
u2, x2, _ = solve(20000,   xn=L, bc0=start_bc, bcn=end_bc, eps=0.0001)

plt.figure()
plt.plot(x4,u4,label="eps=1")
plt.plot(x3,u3, label="eps=0.01")
plt.plot(x2,u2, label="eps=0.0001")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.grid(True)
plt.legend()
plt.show()