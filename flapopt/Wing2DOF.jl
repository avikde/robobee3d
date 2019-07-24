
module Wing2DOF

using LinearAlgebra
include("WingOptimizer.jl")

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
    ρ = 1.225
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
    vaero = [dσ, 0]
    # TODO: confirm and expose design params as argument
    Faero = 1/2 * ρ * cbar * R * (vaero ⋅ vaero) * [CD, CL] * sign(-dσ)

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

    ddq = M \ (-corgrav + τdamp + τaero + τinp)
    # return ddq
    return KSCALE .* [[dσ, dΨ]; ddq]
end

function limits()
    # This is based on observing the OL trajectory
    umin = [-0.15]
    umax = -umin
    xmax = [0.02 * RESCALE, 1.5, Inf, Inf]
    xmin = -xmax
    return umin, umax, xmin, xmax
end

# "Cost function components"
# function eval_f(traj)
#     N = WingOptimizer.getN(traj)
#     [(paero, _, Faero = aero(y); ) for i = 1:N
	
# end

function conEq(x)

end

function conIneq(x)

end

function Jcost(x)

end

# # IPOPT interfacing
# function eval_f(x)
#   return x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
# end

# function eval_g(x, g)
#   g[1] = x[1]   * x[2]   * x[3]   * x[4]
#   g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
# end

# function eval_grad_f(x, grad_f)
#   grad_f[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3])
#   grad_f[2] = x[1] * x[4]
#   grad_f[3] = x[1] * x[4] + 1
#   grad_f[4] = x[1] * (x[1] + x[2] + x[3])
# end

# function eval_jac_g(x, mode, rows, cols, values)
#   if mode == :Structure
#     # Constraint (row) 1
#     rows[1] = 1; cols[1] = 1
#     rows[2] = 1; cols[2] = 2
#     rows[3] = 1; cols[3] = 3
#     rows[4] = 1; cols[4] = 4
#     # Constraint (row) 2
#     rows[5] = 2; cols[5] = 1
#     rows[6] = 2; cols[6] = 2
#     rows[7] = 2; cols[7] = 3
#     rows[8] = 2; cols[8] = 4
#   else
#     # Constraint (row) 1
#     values[1] = x[2]*x[3]*x[4]  # 1,1
#     values[2] = x[1]*x[3]*x[4]  # 1,2
#     values[3] = x[1]*x[2]*x[4]  # 1,3
#     values[4] = x[1]*x[2]*x[3]  # 1,4
#     # Constraint (row) 2
#     values[5] = 2*x[1]  # 2,1
#     values[6] = 2*x[2]  # 2,2
#     values[7] = 2*x[3]  # 2,3
#     values[8] = 2*x[4]  # 2,4
#   end
# end

# function eval_h(x, mode, rows, cols, obj_factor, lambda, values)
#   if mode == :Structure
#     # Symmetric matrix, fill the lower left triangle only
#     idx = 1
#     for row = 1:4
#       for col = 1:row
#         rows[idx] = row
#         cols[idx] = col
#         idx += 1
#       end
#     end
#   else
#     # Again, only lower left triangle
#     # Objective
#     values[1] = obj_factor * (2*x[4])  # 1,1
#     values[2] = obj_factor * (  x[4])  # 2,1
#     values[3] = 0                      # 2,2
#     values[4] = obj_factor * (  x[4])  # 3,1
#     values[5] = 0                      # 3,2
#     values[6] = 0                      # 3,3
#     values[7] = obj_factor * (2*x[1] + x[2] + x[3])  # 4,1
#     values[8] = obj_factor * (  x[1])  # 4,2
#     values[9] = obj_factor * (  x[1])  # 4,3
#     values[10] = 0                     # 4,4

#     # First constraint
#     values[2] += lambda[1] * (x[3] * x[4])  # 2,1
#     values[4] += lambda[1] * (x[2] * x[4])  # 3,1
#     values[5] += lambda[1] * (x[1] * x[4])  # 3,2
#     values[7] += lambda[1] * (x[2] * x[3])  # 4,1
#     values[8] += lambda[1] * (x[1] * x[3])  # 4,2
#     values[9] += lambda[1] * (x[1] * x[2])  # 4,3

#     # Second constraint
#     values[1]  += lambda[2] * 2  # 1,1
#     values[3]  += lambda[2] * 2  # 2,2
#     values[6]  += lambda[2] * 2  # 3,3
#     values[10] += lambda[2] * 2  # 4,4
#   end
# end


# prob = createProblem(n, x_L, x_U, m, g_L, g_U, 8, 10,
#                      eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)

# # Set starting solution
# prob.x = [1.0, 5.0, 5.0, 1.0]

# # Solve
# status = solveProblem(prob)

# println(Ipopt.ApplicationReturnStatus[status])
# println(prob.x)
# println(prob.obj_val)

end
