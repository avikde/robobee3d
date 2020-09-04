import numpy as np
import scipy.linalg

def quadrotorS(Qpos, Qvel, Rdiag=np.ones(4)):
    # For x = (q position, p momentum)
    Md = np.array([100, 100, 100, 3333, 3333, 1000])
    M = np.diag(Md)
    T0 = 9.81e-3 # current thrust
    Z6 = np.zeros((6,6))
    Z3 = np.zeros((3,3))
    dpdotdphi = Md[0] * T0 * np.array([[0,1,0],
            [-1,0,0],
            [0,0,0]])
    dy2dy = np.block([[Z3, dpdotdphi],
        [Z3, Z3]])

    B = np.vstack((
        np.zeros((6,4)), 
        np.array([
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        ))
    A = np.vstack((
        np.hstack((Z6, np.linalg.inv(M))),
        np.hstack((dy2dy, Z6))
    ))
    Q = np.diag(np.hstack((Qpos, Qvel)))
    R = np.diag(Rdiag)
    S = scipy.linalg.solve_continuous_are(A, B, Q, R)
    # aa = A.T @ S + S @ A - S @ B @ np.linalg.inv(R) @ B.T @ S + Q <- should be = 0
    # print(S)
    return S

