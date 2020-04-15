
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once
# using Plots
# Plots.scalefontsizes(0.7) # Only needs to be run once

using BenchmarkTools, Statistics, MAT, Dierckx
using Revise # while developing
import controlutils
cu = controlutils
include("w2d_model.jl")

m = Wing2DOFModel(
	kbo = [30, 0],
	ma = 6,
	ka = 240,
	Amp = deg2rad.([90, 140]))
ny, nu = cu.dims(m)

function getInitialParams(itype=0)
	# robobee scale
	return 75, [3.2^2,  # cbar2[mm^2] (area/R)^2
		2.6666, # τ1 (from 3333 rad/m, [Jafferis (2016)])
		0.55, # mwing[mg] ~=Izz/(mwing*ycp^2). with ycp=8.5, Izz=51.1 [Jafferis (2016)], get
		2.5, # wΨ [mm]
		0, # τ2 quadratic term https://github.com/avikde/robobee3d/pull/92
		54.4, # Aw = 3.2*17 [mm^2] (Jafferis 2016)
		0.0758 # dt
	]
end
uampl, param0 = getInitialParams()
σamax = 0.3 # [mm] constant? for robobee actuators

include("w2d_paramopt.jl") # for opt
N, trajt, traj0, opt, Φ0 = initTraj(m, param0, 1; uampl=uampl)

# # TEST nonlin transmission
# param0[2] = 2
# param0[5] = 2.5 * param0[2]

# "Cost function components" ------------------

skew(a) = [0 -a[3] a[2];
	a[3] 0 -a[1];
	-a[2] a[1] 0]
wrench(paero, Faero) = [Faero; skew(paero) * Faero]
function wrenchy(y, param, flip=false)
	paero, Jaero, Faero = w2daero(m, y, param)
	if flip
		paero[2] = -paero[2]
		Faero[2] = -Faero[2]
	end
	return wrench(paero, Faero)
end

function wrenchAt(inp, param; makeplot=false, Ncycav=5)
	freq, uamplL, dcL, uamplR, dcR, phaseoffs, h2, h3 = inp
	# h2, h3 are 2nd and 3rd harmonics and normalized
	# (1 + h3)sin(wt) + h3 sin(3wt) + h2 sin(2wt)
	Nn = 100
	simdt = 0.02
	yL = createInitialTraj(m, opt, Nn, freq, [1e3, 1e2], param, 0; uampl=uamplL, h3=h3, rawtraj=true, verbose=false, dcoffs=dcL, makeplot=makeplot, simdt=simdt, h2=h2) # ny,N array
	Ntot = size(yL, 2)
	# second harmonic reversed for split cycle yaw
	yR = createInitialTraj(m, opt, Nn, freq, [1e3, 1e2], param, 0; uampl=uamplR, h3=h3, rawtraj=true, verbose=false, dcoffs=dcR, phaseoffs=phaseoffs, makeplot=makeplot, simdt=simdt, h2=-h2)

	# Figure out how many points make up Ncyc cycles
	tNcyc = Ncycav/freq
	Nend = convert(Int64, round(tNcyc/simdt)) # how many to look at
	
	totalWrench = zeros(Nend, 6)
	# loop through the traj, and compute rcop, F and return vector
	for k=1:Nend
		totalWrench[k,:] = wrenchy(yL[:,Ntot-Nend+k], param) + wrenchy(yR[:,Ntot-Nend+k], param, true)
	end
	return mean(totalWrench; dims=1)
end

wrenchNames = ["Fx [mN]", "Fy [mN]", "Fz [mN]", "Rx roll [mN-mm]", "Ry pitch [mN-mm]", "Rz yaw [mN-mm]"]

# wrenchAt([0.16, 75, 0, 75, 0, 0], param0; makeplot=true)

## Single variable ----

# test wrt ampl
# xs = range(40, 80, length=10) # uas
# tws = vcat([wrenchAt([0.16, ua, 0, ua, 0, 0, 0, 0], param0) for ua=xs]...)
# xs = range(-0.5, 0.5, length=10) # ph offs
# tws = vcat([wrenchAt([0.16, 75, 0, 75, 0, p, 0, 0], param0) for p=xs]...)
# xs = range(-0.5, 0.5, length=10) # dc offs
# tws = vcat([wrenchAt([0.16, 75, dc, 75, dc, 0, 0, 0], param0) for dc=xs]...)
xs = range(-0.2, 0.2, length=8) # h2
tws = vcat([wrenchAt([0.16, 75, 0, 75, 0, 0, h2, 0], param0) for h2=xs]...)

pls = [plot(xs, tws[:,c], legend=false, ylabel=wrenchNames[c], lw=2) for c=1:6]
plot(pls...)
gui()

## grids of inputs ------------------

function runInputs(m::Wing2DOFModel, opt, param, freq, uas, phoffs, dcoffs, h2s, h3s)
	i = 0
	Ntotal = length(uas)*length(phoffs)*length(dcoffs)*length(h2s)*length(h3s)
	function runRow(ua, p, dc, h2, h3)
		i+=1
		println(i, "/", Ntotal)
		return [ua p dc h2 h3 wrenchAt([freq, ua+p, dc, ua-p, dc, 0, h2, h3], param0)]
	end
	results = [runRow(ua, p, dc, h2, h3) for ua in uas, p in phoffs, dc in dcoffs, h2 in h2s, h3 in h3s]
	resdict = Dict("results" => results)
	# ^ returns a 2D array result arrays
	matwrite(string("runInputs_", freq, ".zip"), resdict; compress=true)

	return resdict
end
# according to the order in the row
varNames = ["Fact [mN]", "Fact diff [mN]", "DC offs [mN]", "Harmonic 2 [ ]", "Harmonic 3 [ ]"]
for ff = 0.14:0.04:0.22
	# ua, dc
	# runInputs(m, opt, param0, ff, range(40, 80, length=8), [0], range(-20,20, length=8), [0], [0])
	# dc, h2
	# runInputs(m, opt, param0, ff, [75], [0], range(-20,20, length=8), range(-0.2, 0.2, length=8), [0])
	# h2, h3
	# runInputs(m, opt, param0, ff, [75], [0], [0], range(-0.2, 0.2, length=8), range(-0.2, 0.2, length=8))
	# # adiff, dc
	# runInputs(m, opt, param0, ff, [75], range(-20, 20, length=8), range(-20, 20, length=8), [0], [0])
	# ua, adiff
	runInputs(m, opt, param0, ff, range(40, 80, length=8), range(-20, 20, length=8), [0], [0], [0])
end

## Plot against grids of inputs -------------


"ix, iy index into varNames (which is the order of variables in the row)"
function plotInputs(resarg, ix, iy, nvars=5; s=0)
	resdict = typeof(resarg) == String ? matread(resarg) : resarg
	
	# Produced unstructured xi,yi,zi data
	vari = zeros(nvars, 0)
	wri = zeros(6,0)
	for res in resdict["results"]
		# println(res)
		# Append to unstructured data - first the variables
		vari = [vari  res[1:nvars]]
		# then the results
		wri = [wri  res[nvars+1:nvars+6]]
	end
	
	# functions for making splines for plotting
	splineFromUnstructured(xi, yi, zi) = Spline2D(xi, yi, zi; s=s)
	function contourFromUnstructured(xi, yi, zi; xlabel="", ylabel="", title="")
		# If there is no variation at all the spline fails
		if maximum(zi) - minimum(zi) < 5e-2
			return nothing
		end
		spl = splineFromUnstructured(xi, yi, zi)
		function ff(x,y)
			zspl = spl(x,y)
			return zspl >= minimum(zi) && zspl <= maximum(zi) ? zspl : NaN
		end
		# TODO: may have to make these be external
		xpl = (minimum(xi), maximum(xi))
		ypl = (minimum(yi), maximum(yi))
		X = range(xpl..., length=50)
		Y = range(ypl..., length=50)
		pl = contour(X, Y, ff, 
			titlefontsize=10, grid=false, lw=2, c=:bluesreds, 
			xlabel=xlabel, ylabel=ylabel, title=title,
			xlims=xpl, ylims=ypl)
		hline!(pl, [0], ls=:dash, lw=2, color=:black, legend=false)
		return pl
	end

	plsDC = [contourFromUnstructured(vari[ix,:], vari[iy,:], wri[c,:]; xlabel=varNames[ix], ylabel=varNames[iy], title=wrenchNames[c])  for c=1:6]
	filter!(x -> !isnothing(x), plsDC)
	return plsDC
	# 		vcat(
	# 	# [contourFromUnstructured(uai, phoffi, wri[c,:]; xlabel="Fact [mN]", ylabel="phoffs [rad]", title=wrenchNames[c]) 
	# 	# for c=1:6],
	# 	[contourFromUnstructured(uai, dci, wri[c,:]; xlabel="Fact [mN]", ylabel="dcoffs [mN]", title=wrenchNames[c]) 
	# 	for c=1:6]
	# )
end

# pls = plotInputs("runInputs_0.14.zip", 4, 5; s=1000) # h2h3
# pls = plotInputs("runInputs_0.22.zip", 1, 3; s=1000) # ua, dc
# pls = plotInputs("runInputs_0.22.zip", 2, 3; s=1000) # adiff, dc
pls = plotInputs("runInputs_0.14.zip", 1, 2; s=1000) # ua, adiff
# pls = plotInputs("runInputs_0.22.zip", 3, 4; s=1000) # dc, h2
if length(pls) == 2
	plot(pls..., size=(500,250))
else
	plot(pls...)
end
# gui()
savefig("h2h3.png")

## ----

# "Objective to minimize"
# function cu.robj(m::Wing2DOFModel, opt::cu.OptOptions, traj::AbstractArray, param::AbstractArray)::AbstractArray
#     ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    
# 	yk = k -> @view traj[liy[:,k]]
# 	uk = k -> @view traj[liu[:,k]]
    
# 	cbar, τ1, mwing, kΨ, bΨ, τ2, dt = param
	
# 	function wrenchk(k)
# 		paero, _, Faero = w2daero(m, cu.transmission(m, yk(k), param)[1], param)
# 		return wrench(paero, Faero)
# 	end

# 	# Compute the average wrench over the cycle
# 	wrenchAvg = zeros(3)
# 	for k=1:N
# 		wrenchAvg += wrenchk(k)
# 	end
# 	wrenchAvg /= N
#     # avg lift
#     return [wrenchAvg[2] - 100]
# end

# # traj opt ------------------------------------

# εs = [0.05, 0.005, 0.001] # IC, dyn, symm
# prob = cu.ipoptsolve(m, opt, traj0, param0, εs, :traj; print_level=1)
# traj1 = prob.x

# pls = plotTrajs(m, opt, [param0, param0], [traj0, traj1])
# plot(pls...)

# # plot(plot(prob.g), plot(prob.mult_g), size=(900,400))

# mo = nothing#cu.paramoptQPSetup(m, opt, traj0; scaling=false, verbose=false)

# # IPOPT
# εs = [0.05, 0.005, 0.001] # IC, dyn, symm
# prob = cu.ipoptsolve(m, opt, traj0, param0, εs, :traj; print_level=1, nlp_scaling_method="none")
# traj1 = prob.x
# # trajs = [traj0, traj1]
# # params = [param0, param0]

# # # naive param opt
# # εsp = [100.0, 0.005, 100.0] # IC, dyn, symm
# # param1 = cu.ipoptsolve(m, opt, traj0, param0, εs, :param)
# # trajs = [traj0, traj0]
# # params = [param0, param1]

# # with my modification to Ha/Coros g-preferred param opt
# δx = cu.paramδx(m, opt, traj0, param0, prob.mult_x_L, prob.mult_x_U)
# param1 = cu.paramopt(mo, m, opt, traj1, param0, δx, εs, Q; step=1e2)
# # param1 = cu.paramoptJ(m, opt, traj1, param0, εs; step=0.01)
# prob = cu.ipoptsolve(m, opt, traj1, param1, εs, :traj; print_level=1, nlp_scaling_method="none")
# traj2 = prob.x

# trajs = [traj0, traj1, traj2]
# params = [param0, param0, param1]

# # # Custom solver ---
# # wkt = cu.OptWorkspace(cu.Ntraj(m, opt, N), (N+2)*ny)
# # @time traj1 = cu.csSolve!(wkt, m, opt, traj0, param0, :traj; Ninner=30, μs=[1e5])
# # wkp = cu.OptWorkspace(length(param0), (N+2)*ny)
# # @time param1 = cu.csSolve!(wkp, m, opt, traj1, param0, :param; Ninner=30, μs=[1e3])
# # @time traj2 = cu.csSolve!(wkt, m, opt, traj1, param1, :traj; Ninner=30, μs=[1e5])
# # trajs = [traj0, traj1, traj2]
# # params = [param0, param0, param1]
# # # trajs, params, wkt = cu.csAlternateSolve(m, opt, traj0, params0, 1; μst=[1e6], Ninnert=30, μsp=[1e-2,1e-2], Ninnerp=2)

# # # pl2 = plotParams(m, opt, trajs[:,end], (params[:,i] for i = 1:size(params,2))...; μ=1e-1)
# # # display(params)

# # paramopt -------------------------------------


# # visualize --------------------------------------------

# println("Objectives: ", [(_ro = cu.robj(m, opt, tt, param0); _ro ⋅ _ro) for tt in trajs])
# println("Params: ", params)

# # animateTrajs(m, opt, params, trajs)
# pl1 = plotTrajs(m, opt, params, trajs)
# pl2 = cu.visualizeConstraintViolations(m, opt, params, trajs)

# l = @layout [grid(2,2) a]
# plot(pl1..., pl2, layout=l, size=(900,400))

# # plot(pl1...)
# # gui()
