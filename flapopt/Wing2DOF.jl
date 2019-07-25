
module Wing2DOF

using LinearAlgebra, StaticArrays, DifferentialEquations
include("DirTranForm.jl")

const RESCALE = 30.0
const KSCALE = @SVector [RESCALE, 1, RESCALE, 1]

const ny = 4
const nu = 1

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
    σ, Ψ, dσ, dΨ = (@SVector [T, 1, T, 1]) .* y
    cΨ = cos(Ψ)
    sΨ = sin(Ψ)
    α = π / 2 - Ψ # AoA

    # aero force
    wing1 = @SVector [σ, 0]
    RΨ = @SMatrix [cΨ -sΨ; sΨ cΨ]
    paero = wing1 + RΨ * @SVector [0, -cbar]
    Jaero = @SMatrix [1 cbar * cΨ; 0 cbar * sΨ]
    Caero = @SVector [(CDmax + CD0)/2 - (CDmax - CD0)/2 * cos(2 * α), CLmax * sin(2 * α)]
    vaero = @SVector [dσ, 0]
    Faero = 1/2 * ρ * cbar * R * (vaero ⋅ vaero) * Caero * sign(-dσ)

    return paero, Jaero, Faero
end

"Continuous dynamics second order model"
function dydt(yin::Vector, u::Vector, _params::Vector)
    # unpack
    cbar, T = _params
    y = (1 ./ KSCALE) .* yin
    σ, Ψ, dσ, dΨ = (@SVector [T, 1, T, 1]) .* y
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
    M = @SMatrix [mspar+mwing   cbar*mwing*cΨ; cbar*mwing*cΨ   Iwing+cbar^2*mwing]
    corgrav = @SVector [ka*σ - cbar*mwing*sΨ*dΨ^2, khinge*Ψ]
    # non-lagrangian terms
    τdamp = @SVector [0, -bΨ * dΨ]
    _, Jaero, Faero = aero(y, u, _params)
    τaero = Jaero' * Faero
    # input
    τinp = @SVector [u[1]/T, 0]

    ddq = inv(M) * (-corgrav + τdamp + τaero + τinp)
    # return ddq
    return KSCALE .* [dσ, dΨ, ddq[1], ddq[2]]
end

function limits()
    # This is based on observing the OL trajectory
    umin = @SVector [-0.15]
    umax = -umin
    xmax = @SVector [0.02 * RESCALE, 1.5, Inf, Inf]
    xmin = -xmax
    return umin, umax, xmin, xmax
end

# Create an initial traj --------------

function createInitialTraj(freq::Real, posGains::Vector, params0::Vector)
    # Create a traj
    σmax = Wing2DOF.limits()[end][1]
    function strokePosController(y, t)
        σdes = 0.75 * σmax * sin(freq * 2 * π * t)
        return posGains[1] * (σdes - y[1]) - posGains[2] * y[3]
    end
    strokePosControlVF(y, p, t) = Wing2DOF.dydt(y, [strokePosController(y, t)], params0)
    # OL traj1
    teval = collect(0:1e-4:0.1)
    prob = ODEProblem(strokePosControlVF, zeros(4), (teval[1], teval[end]))
    sol = solve(prob, saveat=teval)

    # σt = plot(sol, vars=1, linewidth=1, ylabel="act disp [m]")
    # Ψt = plot(sol, vars=2, linewidth=1, ylabel="hinge ang [r]")
    # plot(σt, Ψt, layout=(2,1))
    # gui()

    olRange = 171:3:238
    trajt = sol.t[olRange]
    δt = trajt[2] - trajt[1]
    N::Int = length(trajt) - 1
    olTrajaa = sol.u[olRange] # 23-element Array{Array{Float64,1},1} (array of arrays)
    olTraju = [strokePosController(olTrajaa[i], trajt[i]) for i in 1:N] # get u1,...,uN
    traj0 = [vcat(olTrajaa...); olTraju; δt] # dirtran form {x1,..,x(N+1),u1,...,u(N),δt}

    return trajt, traj0
end

# "Cost function components" ------------------

"Objective, min"
function eval_f(traj, params)
    N = DirTranForm.getN(traj, ny, nu)
    Favg = @SVector zeros(2)
    ly, lu = DirTranForm.linind(traj, ny, nu)
    for k = 1:N
        vy = @view ly[:,k]
        vu = @view lu[:,k]
        paero, _, Faero = aero(traj[vy], traj[vu], params)
        Favg += Faero
    end
    # max avg lift
	return -Favg[2]
end

"Inequality and equality constraints"
function eval_g!(traj, g, params)
    N = DirTranForm.getN(traj, ny, nu)
    ly, lu = DirTranForm.linind(traj, ny, nu)
    δt = traj[end]
    for k = 1:N
        vy = @view ly[:,k]
        vy2 = @view ly[:,k+1]
        vu = @view lu[:,k]
        g[vy] = traj[vy2] - (traj[vy] + δt * dydt(traj[vy], traj[vu], params))
    end
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
