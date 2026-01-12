#h = 3
import numpy as np
from scipy.sparse import diags
import time

t0 = time.time()

def K_matrix(h):
    k_11 = 1/h + h/3
    k_12 = -1/h + h/6
    return k_11, k_12

def A_upper_algo_1(h_n):
    n = h_n.shape[0]+1
    A = np.zeros([n,n])
    for i in range(n-1):
        h = h_n[i]
        k_11, k_12 = K_matrix(h)
        A[i,i] += k_11
        A[i,i+1] = k_12
        A[i+1,i+1] = k_11
    return A

def A_upper_algo_2(A,c,d):
    n = A_upper.shape[0]
    b = np.zeros(n)
    b[0] = c
    b[1] = b[1]-A[0,1]*c

    A[0,0] = 1
    A[0,1] = 0

    b[n-1] = d
    b[n-2] = b[n-2]-A[n-2,n-1]*d
    A[n-1,n-1] = 1
    A[n-2,n-1] = 0
    return A,b

c = 1
d = np.exp(2)
h_cons = np.zeros(100)+1

A_upper = A_upper_algo_1(h_cons)
A_upper,b = A_upper_algo_2(A_upper,c,d)

A_sparse = diags(
    diagonals=[A_upper.diagonal(1), A_upper.diagonal(), A_upper.diagonal(1)],
    offsets=[-1, 0, 1],
    format="csr"
)
A_full = A_sparse.toarray()

np.linalg.solve(A_full,b)
t1 = time.time()
print(t1-t0)