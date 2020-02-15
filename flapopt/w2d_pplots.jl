using Plots

"Param space plots"
function plotParams(m::Wing2DOFModel, opt::cu.OptOptions, rr; compareTo=nothing)
    # # First plot the param landscape
    # pranges = [
    #     0:0.25:6.0, # cbars
    #     10.0:1.0:50, # Ts
    #     0.1:0.1:3.0, # mwings
    #     0.1:1.0:20.0, # kΨs
    #     0.1:1.0:20.0 # bΨs
    # ]
    # labels = [
    #     "chord",
    #     "T",
    #     "mwing",
    #     "hinge k",
    #     "hinge b"
    # ]

    # # different param vectors passed in
    # params = hcat(args...) # Np x Nsteps
	# param0 = args[end] # for the slices use the last param
	
    function plotSlice(i1, i2, xlabel, xpl, ylabel, ypl)
        function f(p1, p2)
            x = copy(rr["x"])
            x[i1] = p1
            x[i2] = p2
            return rr["eval_f"](x)
        end
        X = range(xpl..., length=50)
        Y = range(ypl..., length=50)
        pl = contour(X, Y, f, 
			titlefontsize=10, grid=false, lw=2, c=:bluesreds, 
			xlabel=xlabel, ylabel=ylabel, #title=title,
			xlims=xpl, ylims=ypl)
		# Now plot the locations of different params
		plparam(p; kwargs...) = plot!(pl, [p[i1]], [p[i2]], legend=false, markershape=:auto; kwargs...)
		plparam(rr["param"]; markercolor=:red)
		if !isnothing(compareTo)
			plparam(compareTo; markercolor=:green)
        end

        # Test superimpose constraint
        ftest(Aw, dt) = Aw - 80
        contour!(pl, X, Y, ftest, fill=false)
        return pl
    end
	
	return [
		plotSlice(6, 7, "Aw", (50, 100), "dt", (0.05, 0.1)),
		plotSlice(2, 7, "T1", (1.6, 3.2), "dt", (0.05, 0.1))
		]
end

# function plotParamImprovement(m::Wing2DOFModel, opt::cu.OptOptions, params, trajs, paramObj::Function)
#     ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, trajs[1])

#     # The param space plots
#     pls = plotParams(m, opt, trajs[1], paramObj, params...)
#     # Traj plots
#     σt, Ψt, ut, liftt, dragt = plotTrajs(m, opt, params, trajs)

#     return pls..., σt, Ψt, ut, liftt, dragt
# end

function debugConvexity(m, opt, POPTS, ret1)
    # Turn down wΔy = 1e2 to make weight be just for LSE https://github.com/avikde/robobee3d/pull/140#issuecomment-583749876
    wΔy0 = copy(POPTS.R[4])
    dtlimU0 = copy(POPTS.plimsU[end])
    println("HI1 ", POPTS.R[4])
    # POPTS.R[4] = 1e2
    
    # Qdt = 0 -> f = 151, FD = 60, J=83.  but making f >= 0.18 gives 
    ret2 = @time opt1(m, ret1["traj"], ret1["param"], 1, 180; Φ=90, Qdt=0.0, tol=5e-2)

    POPTS.plimsU[end] = 1/(N*0.18)

    ret3 = @time opt1(m, ret1["traj"], ret1["param"], 1, 180; Φ=90, Qdt=0.0, tol=5e-2)

    # replace
    # POPTS.R[4] = wΔy0
    POPTS.plimsU[end] = dtlimU0

    J(ret) = ret["eval_f"](ret["x"])
    println("HIHI ", J(ret2), "->", J(ret3))

    ##
    # pls = plotParams(m, opt, ret2; compareTo=[16.365 2.462 0.85 4.084 4.924 65.457 0.069])
    # pls = plotParams(m, opt, ret2; compareTo=[22.109 2.462 1.149 4.415 4.924 88.436 0.094])
    # plot(pls...)
    # gui()
end
