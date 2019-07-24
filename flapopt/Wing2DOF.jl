
module Wing2DOF

using LinearAlgebra

# nx = 4
# nu = 1
# y0 = np.zeros(nx)

const RESCALE = 30.0
const KSCALE = [RESCALE, 1, RESCALE, 1]

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

"Continuous dynamics second order model"
function dydt(yin::Vector, u::Vector, _params::Vector)
    # unpack
    cbar, T = _params
    y = (1 ./ KSCALE) .* yin
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
    τdamp = [0, -bΨ * dΨ]
    _, Jaero, Faero = aero(y, u, _params)
    τaero = Jaero' * Faero
    # input
    τinp = [u[1]/T, 0]

    ddq = inv(M) * (-corgrav + τdamp + τaero + τinp)
    # return ddq
    return KSCALE .* [[dσ, dΨ]; ddq]
end

# @property
# def limits(self):
#     # This is based on observing the OL trajectory
#     umin = np.array([-0.15])
#     umax = -umin
#     xmin = np.array([-0.02 * self.rescale, -1.2, -np.inf, -np.inf])
#     xmax = -xmin
#     return umin, umax, xmin, xmax

end
