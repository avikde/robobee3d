"Param space plots"
include("w2d_model.jl")

function plotParams(m::Wing2DOFModel, opt::cu.OptOptions, traj::Vector, paramObj::Function, args...; μ::Float64=1e-1, Vclip=50000)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    # First plot the param landscape
    pranges = [
        0:0.25:6.0, # cbars
        10.0:1.0:50, # Ts
        0.1:0.1:3.0, # mwings
        0.1:1.0:20.0, # kΨs
        0.1:1.0:20.0 # bΨs
    ]
    labels = [
        "chord",
        "T",
        "mwing",
        "hinge k",
        "hinge b"
    ]

    # different param vectors passed in
    params = hcat(args...) # Np x Nsteps
    param0 = args[end] # for the slices use the last param

    # # Old: f defined here
    # Ng = opt.boundaryConstraint == cu.SYMMETRIC ? (N+2)*ny : (N+1)*ny
    # g = zeros(Ng)#cu.gbounds(m, opt, traj)[1]
    # f(p1, p2) = begin
    #     cu.gvalues!(g, m, opt, traj, [p1,p2], traj[1:4])
    #     cu.Jobj(m, opt, traj, [p1,p2]) + μ/2 * g' * g
    # end

    function plotSlice(i1, i2)
        function f(p1, p2)
            parg = copy(param0)
            parg[i1] = p1
            parg[i2] = p2
            V = paramObj(parg)
            return V > Vclip ? NaN : V
        end
        pl = contour(pranges[i1], pranges[i2], f, fill=true, seriescolor=cgrad(:bluesreds), xlabel=labels[i1], ylabel=labels[i2])
        # Now plot the path taken
        plot!(pl, params[i1,:], params[i2,:], marker=:auto, legend=false)
        # just in case
        xlims!(pl, (pranges[i1][1], pranges[i1][end]))
        ylims!(pl, (pranges[i2][1], pranges[i2][end]))
        return pl
    end
    
    return (plotSlice(1, 2), # cbar, T slice
        plotSlice(1, 3), # cbar, mwing slice
        plotSlice(2, 3), # T, mwing slice
        plotSlice(4, 5) # T, khinge slice
    )
end

function plotParamImprovement(m::Wing2DOFModel, opt::cu.OptOptions, params, trajs, paramObj::Function)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, trajs[1])

    # The param space plots
    pls = plotParams(m, opt, trajs[1], paramObj, params...)
    # Traj plots
    σt, Ψt, ut, liftt, dragt = plotTrajs(m, opt, params, trajs)

    return pls..., σt, Ψt, ut, liftt, dragt
end
