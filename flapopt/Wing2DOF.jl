using LinearAlgebra
using ForwardDiff


# nx = 4
# nu = 1
# y0 = np.zeros(nx)

RESCALE = 30

"""
Returns aero force
"""
function aero(y::Vector, u::Vector, _params::Vector)
    CLmax = 1.8
    CDmax = 3.4
    CD0 = 0.4
    rho = 1.225
    R = 15e-3
    
    # unpack
    cbar, T = _params
    σ, Ψ, dσ, dΨ = [T, 1, T, 1] .* y
    cΨ = cos(Ψ)
    sΨ = sin(Ψ)
    α = π / 2 - Ψ # AoA

    # aero force
    wing1 = [σ, 0]
    paero = wing1 + [cΨ -sΨ; sΨ cΨ] * [0, -cbar]
    Jaero = [1 cbar * cΨ; 0 cbar * sΨ]
    CL = CLmax * sin(2 * α)
    CD = (CDmax + CD0)/2 - (CDmax - CD0)/2 * cos(2 * α)
    vaero = [dΨ, 0]
    # TODO: confirm and expose design params as argument
    Faero = 1/2 * rho * cbar * R * dot(vaero, vaero) * [CD, CL] * sign(-dΨ)

    return paero, Jaero, Faero
end

# "Continuous dynamics"
function dydt(yin::Vector, u::Vector, _params::Vector)
    y = [1/RESCALE, 1, 1/RESCALE, 1] .* yin
    # unpack
    cbar, T = _params
    σ, Ψ, dσ, dΨ = [T, 1, T, 1] .* y
    # NOTE: for optimizing transmission ratio
    # Thinking of y = (sigma_actuator, psi, dsigma_actuator, dpsi)
    # u = (tau_actuator)
    # sigma = sigma_actuator * T; tau = tau_actuator / T
    cΨ = cos(Ψ)
    sΨ = sin(Ψ)

    # params
    mspar = 0
    ka = 0
    khinge = 1e-3
    mwing = 5e-6
    Iwing = 1e-9#mwing * cbar**2
    bΨ = 5e-7

    # inertial terms
    M = [mspar+mwing   cbar*mwing*cΨ; cbar*mwing*cΨ   Iwing+cbar^2*mwing]
    corgrav = [ka*σ - cbar*mwing*sΨ*dΨ^2, khinge*Ψ]
    # non-lagrangian terms
    taudamp = [0, -bΨ * dΨ]
    _, Jaero, Faero = aero(y, u, _params)
    tauaero = Jaero' * Faero
    # input
    tauinp = [u[1]/T, 0]

    ddq = inv(M) * (-corgrav + taudamp + tauaero + tauinp)

    return RESCALE .* ddq
end

println("hi")
y0 = [0.1,0.1,1,0]
u0 = [0.]
params0 = [0.05,1.0]
paeroFun(q::Vector) = aero([q;[0,0]], u0, params0)[1]
println("paero ", paeroFun(y0[1:2]))

# Try to use ForwardDiff
JaeroFun = y -> ForwardDiff.jacobian(paeroFun, y)
# println(JaeroFun(y0[1:2]))
println(aero(y0, u0, params0)[2])


println(dydt(y0, u0, params0))

# @property
# def limits(self):
#     # This is based on observing the OL trajectory
#     umin = np.array([-0.15])
#     umax = -umin
#     xmin = np.array([-0.02 * self.rescale, -1.2, -np.inf, -np.inf])
#     xmax = -xmin
#     return umin, umax, xmin, xmax