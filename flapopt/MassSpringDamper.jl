
using DifferentialEquations
using Plots; gr()

import controlutils
cu = controlutils

struct MassSpringDamperModel <: controlutils.Model end

# gear ratio
const G = 1.0
const mb = 1.0

function cu.dims(m::MassSpringDamperModel)::Tuple{Int, Int}
    return 2, 1
end

function cu.pdims(m::MassSpringDamperModel)::Int
    return 1
end

function cu.plimits(m::MassSpringDamperModel)
    return [0.1, 1.0]
end

function cu.limits(m::MassSpringDamperModel)::Tuple{Vector, Vector, Vector, Vector}
    # This is based on observing the OL trajectory. See note on units above.
    umax = [75.0] # [mN]
    umin = -umax
    xmax = [10,1000] # [mm, rad, mm/ms, rad/ms]
    xmin = -xmax
    return umin, umax, xmin, xmax
end

function cu.limitsTimestep(m::MassSpringDamperModel)::Tuple{Float64, Float64}
	return 0.01, 0.07
end

"Continuous dynamics second order model"
function cu.dydt(model::MassSpringDamperModel, y::AbstractArray, u::AbstractArray, _params::Vector)::AbstractArray
    # unpack
    kσ = _params[1]
    σ, σ̇ = [1/G, 1/G] .* y # [mm, mm/ms]
    # NOTE: for optimizing transmission ratio
    # Thinking of y = (sigma_actuator, psi, dsigma_actuator, dpsi)
    # u = (tau_actuator)
	# sigma = sigma_actuator * T; tau = tau_actuator / T
	
    bσ = 1 # [mN/(mm/ms)]
	
    ddq = 1.0/mb * (-kσ*σ - bσ*σ + G*u[1])
    # return ddq
    return [y[2], G * ddq]
end


"""
freq [kHz]; posGains [mN/mm, mN/(mm-ms)]; [mm, 1]
Example: trajt, traj0 = Wing2DOF.createInitialTraj(0.15, [1e3, 1e2], params0)
"""
function createInitialTraj(m::MassSpringDamperModel, opt::cu.OptOptions, N::Int, freq::Real, posGains::Vector, params::Vector)
    # Create a traj
    σmax = cu.limits(m)[end][1]
    function controller(y, t)
        σdes = 0.9 * σmax * sin(freq * 2 * π * t)
        return posGains[1] * (σdes - y[1]) - posGains[2] * y[2]
    end
    vf(y, p, t) = cu.dydt(m, y, [controller(y, t)], params)
    # OL traj1
    simdt = 0.1 # [ms]
    teval = collect(0:simdt:100) # [ms]
    prob = ODEProblem(vf, [1.,0.], (teval[1], teval[end]))
    sol = solve(prob, saveat=teval)

    # # Animate whole traj
    # Nt = length(sol.t)
    # @gif for k = 1:3:Nt
    #     yk = sol.u[k]
    #     uk = [strokePosController(yk, sol.t[k])]
    #     drawFrame(m, yk, uk, params)
    # end
    # Plot
    σt = plot(sol, vars=1, ylabel="act pos [m/s]")
    plot(σt)
    gui()

    # expectedInterval = opt.boundaryConstraint == cu.SYMMETRIC ? 1/(2*freq) : 1/freq # [ms]
    # expectedPts = expectedInterval / simdt

    # starti = 170
    # olRange = starti:2:(starti + 2*N)
    # trajt = sol.t[olRange]
    # δt = trajt[2] - trajt[1]
    # olTrajaa = sol.u[olRange] # 23-element Array{Array{Float64,1},1} (array of arrays)
    # olTraju = [strokePosController(olTrajaa[i], trajt[i]) for i in 1:N] # get u1,...,uN
    # traj0 = [vcat(olTrajaa...); olTraju] # dirtran form {x1,..,x(N+1),u1,...,u(N),δt}
    # if opt.vart
    #     traj0 = [traj0; δt]
    # else
    #     println("Initial traj δt=", δt, ", opt.fixedδt=", opt.fixedδt)
    # end

    # return trajt .- trajt[1], traj0
end
