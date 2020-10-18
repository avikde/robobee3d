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

def updateConstraint(A, dt, T0, s0s, Btaus):
    # Left third
    AxidxT0dt = []
    n2 = 2*ny + 3 # nnz in each block col on the left
    for k in range(N-2):
        AxidxT0dt += [n2*k + i for i in [8,11,14]]

    # Middle third
    n1 = (2*N-1)*ny + (N-2)*3 # All the nnz in the left third
    n2 = 3*ny # nnz in each of the first N-1 block cols in the middle third
    Axidxdt = []
    for k in range(N):
        if k < N-1:
            Axidxdt += [n1 + n2*k + i for i in [0,3,6,9,12,15]]
        else:
            Axidxdt += [n1 + n2*k + i for i in [0,2,4,6,8,10]]
    
    # Right third
    n1 += 3*ny*(N-1) + 2*ny # all nnz in the left and middle third
    n2 = 9 # nnz in each B0
    Axidxs0 = []
    AxidxBtau = []
    for k in range(N):
        Axidxs0 += [n1 + n2*k + i for i in range(3)]
        AxidxBtau += [n1 + n2*k + 3 + i for i in range(6)]

    # Last check
    assert A.nnz == n1 + n2*N

    # Update
    A.data[AxidxT0dt] = dt*T0
    A.data[Axidxdt] = dt
    A.data[Axidxs0] = dt*np.hstack((s0s))
    A.data[AxidxBtau] = dt*np.hstack([np.ravel(Btau,order='F') for Btau in Btaus])

    # print(A[:,36:42].toarray())

def testDynamicsConstraint():
    ys = np.zeros((N,ny))
    dys = np.zeros((N,ny))
    us = np.random.rand(N,nu)

    y0 = np.random.rand(ny)
    dy0 = np.random.rand(ny)
    T0 = 1
    dt = 1
    
    # Basic constituents of dynamics A0, B0 (only sparsity matters)
    B0 = np.block([
        [np.ones((3,1)),np.zeros((3,2))],
        [np.zeros((3,1)),np.ones((3,2))]
    ])
    A0 = np.zeros((6, 6))
    A0[:3,3:] = np.eye(3)
    c0 = np.ones(ny) # FIXME:

    y1 = y0 + dt * dy0
    yy = np.copy(y1)
    dyy = np.copy(dy0)
    for k in range(N):
        yy1 = yy + dt * dyy
        dyy1 = dyy + (A0 @ yy + B0 @ us[k,:] + c0)

        dys[k,:] = dyy1
        ys[k,:] = yy1

        yy = yy1
        dyy = dyy1

    # stack
    x = np.hstack((np.ravel(ys), np.ravel(dys), np.ravel(us)))
    print(A @ x)

if __name__ == "__main__":
    T0 = 0.5
    dt = 0.5
    s0s = [[0.1,0.1,0.9] for i in range(N)]
    Btaus = [np.full((3,2),1.123) for i in range(N)]

    A = initConstraint()
    updateConstraint(A, dt, T0, s0s, Btaus)
    testDynamicsConstraint()
