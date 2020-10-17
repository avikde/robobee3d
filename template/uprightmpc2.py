import osqp
import numpy as np
import scipy.sparse as sp

ny = 6
nu = 3
N = 3
nx = N * (2*ny + nu)
nc = 2*N*ny

def initConstraint():
    # these will be updated
    T0 = 1
    dt = 1

    A = np.zeros((nc, nx))
    # l = np.zeros(nc)
    # u = np.zeros(nc)

    # cols of A are broken up like this (partition nx)
    n1 = N*ny
    n2 = 2*N*ny
    # rows of A broken up:
    nc1 = N*ny

    # Basic constituents of dynamics A0, B0 (only sparsity matters)
    B0 = np.block([
        [np.ones((3,1)),np.zeros((3,2))],
        [np.zeros((3,1)),np.ones((3,2))]
    ])
    A0 = np.zeros((6, 6))
    A0[:3,3:] = np.eye(3)
    
    for k in range(N):
        # ykp1 = yk + dt*dyk equations
        A[k*ny:(k+1)*ny, k*ny:(k+1)*ny] = -np.eye(ny)
        A[k*ny:(k+1)*ny, n1 + k*ny:n1 + (k+1)*ny] = dt * np.eye(ny)
        if k>0:
            A[k*ny:(k+1)*ny, (k-1)*ny:(k)*ny] = np.eye(ny)
        
        # dykp1 equation
        A[nc1 + k*ny:nc1 + (k+1)*ny, n1 + k*ny:n1 + (k+1)*ny] = -np.eye(ny)
        A[nc1 + k*ny:nc1 + (k+1)*ny, n2 + k*nu:n2 + (k+1)*nu] = B0
        if k>0:
            A[nc1 + k*ny:nc1 + (k+1)*ny, n1 + (k-1)*ny:n1 + (k)*ny] = np.eye(ny)
        if k>1:
            A[nc1 + k*ny:nc1 + (k+1)*ny, (k-2)*ny:(k-1)*ny] = A0
    
    return sp.csc_matrix(A)

if __name__ == "__main__":
    A = initConstraint()
    print(A)
