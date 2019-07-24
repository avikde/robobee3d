using LinearAlgebra
using ForwardDiff


# nx = 4
# nu = 1
# y0 = np.zeros(nx)

# rescale = 1.0

"""
Returns aero force
"""
function aero(y::Vector, u::Vector, _params::Vector)
    cbar = _params[1]
    T = _params[2]
    CLmax = 1.8
    CDmax = 3.4
    CD0 = 0.4
    rho = 1.225
    R = 15e-3
    
    sigma = y[1] * T
    psi = y[2]
    dsigma = y[3] * T
    dpsi = y[4]
    cpsi = cos(psi)
    spsi = sin(psi)
    alpha = pi / 2 - psi

    # aero force
    wing1 = [sigma, 0]
    paero = wing1 + [cpsi -spsi; spsi cpsi] * [0, -cbar]
    Jaero = [1 cbar * cpsi; 0 cbar * spsi]
    CL = CLmax * sin(2 * alpha)
    CD = (CDmax + CD0)/2 - (CDmax - CD0)/2 * cos(2 * alpha)
    vaero = [dsigma, 0]
    # TODO: confirm and expose design params as argument
    Faero = 1/2 * rho * cbar * R * dot(vaero, vaero) * [CD, CL] * sign(-dsigma)

    return paero, Jaero, Faero
end

println("hi")
y0 = [0.1,0.1,1,0]
u0 = [0.]
params0 = [0.05,1]
# paero, Jaero, Faero = aero(y, u, params)
paeroFun(q::Vector) = aero([q;[0,0]], u0, params0)[1]
println("paero ", paeroFun(y0))

# Try to use ForwardDiff
JaeroFun = y -> ForwardDiff.jacobian(paeroFun, y)
println(JaeroFun(y0[1:2]))
println(aero(y0, u0, params0)[2])

# function dydt(self, yin, u, _params=[], dydt)
#     ''' 
#     See mma file flapping wing traj
#     '''
#     cbar = _params[0]
#     T = _params[1]
#     Kscale = np.diag([self.rescale, 1, self.rescale, 1])
#     y = np.linalg.inv(Kscale) @ yin
#     # NOTE: for optimizing transmission ratio
#     # Thinking of y = (sigma_actuator, psi, dsigma_actuator, dpsi)
#     # u = (tau_actuator)
#     # sigma = sigma_actuator * T; tau = tau_actuator / T
#     sigma = y[0] * T
#     psi = y[1]
#     dsigma = y[2] * T
#     dpsi = y[3]
#     cpsi = np.cos(psi)
#     spsi = np.sin(psi)

#     # params
#     mspar = 0
#     ka = 0
#     khinge = 1e-3
#     mwing = 5e-6
#     Iwing = 1e-9#mwing * cbar**2
#     bpsi = 5e-7

#     # inertial terms
#     M = np.array([[mspar + mwing, cbar * mwing * cpsi], [cbar * mwing * cpsi, Iwing + cbar**2 * mwing]])
#     corgrav = np.array([ka * sigma - cbar * mwing * spsi * dpsi**2, khinge * psi])
#     # non-lagrangian terms
#     taudamp = np.array([0, -bpsi * dpsi])
#     _, Jaero, Faero = self.aero(y, u, params)
#     tauaero = Jaero.T @ Faero
#     # input
#     tauinp = np.array([u[0] / T, 0])

#     ddq = np.linalg.inv(M) @ (-corgrav + taudamp + tauaero + tauinp)

#     dydt = Kscale @ np.array([dsigma, dpsi, ddq[0], ddq[1]])
# end

# @property
# def limits(self):
#     # This is based on observing the OL trajectory
#     umin = np.array([-0.15])
#     umax = -umin
#     xmin = np.array([-0.02 * self.rescale, -1.2, -np.inf, -np.inf])
#     xmax = -xmin
#     return umin, umax, xmin, xmax