
module Wing2DOF

using LinearAlgebra, StaticArrays, DifferentialEquations
using Plots; gr()
include("DirTranForm.jl")

const RESCALE = 1.0 # For realistic params
const KSCALE = @SVector [RESCALE, 1, RESCALE, 1]

const ny = 4
const nu = 1
const R = 20e-3

#=
2DOF wing model config:
- q[1] = actuator displacement (m)
- q[2] = hinge angle
- y = (q, dq)
- u = [τ], where τ = actuator force

Reasonable limits:
- MRL Bee: umax = 75mN, y[1]max = 150 μm, mass = 25 mg
- MRL HAMR: umax = 270mN, y[1]max = 123 μm, mass = 100 mg [Doshi et. al. (2015)]

Params:
- cbar = wing chord
- T = transmission ratio. i.e. wing spar displacement σ = T q[1].
- R = wing length. This is irrelevant for the 2DOF model, but choosing one for now as 20mm. This is used to (a) calculate the aero force, and (b) in the plotting to show the "stroke angle" as a more intuitive version of σ.

NOTE:
- The "T" above is unitless. You can intuitively think of a wing "stroke angle" ~= T q[1] / (R / 2) (using σ as the arc length and "R/2" as the radius). This is distinct from the Toutput (in rad/m), and they are related as T ~= Toutput ⋅ (R / 2).
- Reasonable values: Toutput = 2666 rad/m in the current two-wing vehicle; 2150 rad/m in X-Wing; 3333 rad/m in Kevin Ma's older vehicles. With the R above, this suggests T ~= 20-30.
=#

"Returns aero force"
function aero(y::Vector, u::Vector, _params::Vector)
    CLmax = 1.8
    CDmax = 3.4
    CD0 = 0.4
    ρ = 1.225
    
    # unpack
    cbar, T = _params
    σ, Ψ, dσ, dΨ = (@SVector [T, 1, T, 1]) .* y
    cΨ = cos(Ψ)
    sΨ = sin(Ψ)

    # CoP kinematics
    wing1 = @SVector [σ, 0]
    RΨ = @SMatrix [cΨ -sΨ; sΨ cΨ] # Ψ > 0 => hinge looks like \; Ψ < 0 => hinge looks like /
    paero = wing1 + RΨ * @SVector [0, -cbar]
    # FIXME: add moving CoP
    Jaero = @SMatrix [1 cbar * cΨ; 0 cbar * sΨ]
    
    # Aero force
    #=
    Use this in Mathematica to debug the signs.
    Manipulate[
    \[Alpha] = -\[Pi]/2 + \[Psi];
    Graphics@
    Arrow[{{0, 0},
        {((CDmax + CD0)/2 - (CDmax - CD0)/2*Cos[2 \[Alpha]]), CLmax Sin[2 \[Alpha]]} Sign[-d\[Sigma]]}],
    {d\[Sigma], -1, 1}, {\[Psi], -\[Pi]/2, \[Pi]/2}]
    =#
    α = π/2 - Ψ # AoA
    Caero = @SVector [((CDmax + CD0)/2 - (CDmax - CD0)/2 * cos(2α)), CLmax * sin(2α)]
    Faero = 1/2 * ρ * cbar * R * dσ^2 * Caero * sign(-dσ)

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
    mwing = 1e-6
    Iwing = mwing * cbar^2
    kσ = 1e-10
    kΨ = 0
    bΨ = 1e-7

    # inertial terms
    M = @SMatrix [mspar+mwing   cbar*mwing*cΨ; cbar*mwing*cΨ   Iwing+cbar^2*mwing]
    corgrav = @SVector [kσ*σ - cbar*mwing*sΨ*dΨ^2, kΨ*Ψ]
    # non-lagrangian terms
    τdamp = @SVector [0, -bΨ * dΨ]
    _, Jaero, Faero = aero(y, u, _params)
    τaero = Jaero' * Faero
    # input
    τinp = @SVector [u[1]/T, 0]

    ddq = inv(M) * (-corgrav + τdamp + τaero + τinp)
    # return ddq
    return KSCALE .* [y[3], y[4], ddq[1], ddq[2]]
end

function limits()
    # This is based on observing the OL trajectory
    umax = @SVector [75e-3]
    umin = -umax
    xmax = @SVector [150e-6 * RESCALE, 1.5, Inf, Inf]
    xmin = -xmax
    return umin, umax, xmin, xmax
end

# Create an initial traj --------------

function createInitialTraj(freq::Real, posGains::Vector, params0::Vector)
    # Create a traj
    σmax = Wing2DOF.limits()[end][1]
    function strokePosController(y, t)
        σdes = 0.9 * σmax * sin(freq * 2 * π * t)
        return posGains[1] * (σdes - y[1]) - posGains[2] * y[3]
    end
    strokePosControlVF(y, p, t) = Wing2DOF.dydt(y, [strokePosController(y, t)], params0)
    # OL traj1
    teval = collect(0:1e-4:0.1)
    prob = ODEProblem(strokePosControlVF, zeros(4), (teval[1], teval[end]))
    sol = solve(prob, saveat=teval)

    σt = plot(sol, vars=1, linewidth=1, ylabel="act disp [m]")
    Ψt = plot(sol, vars=2, linewidth=1, ylabel="hinge ang [r]")
    plot(σt, Ψt, layout=(2,1))
    gui()

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
