
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
using SparseArrays, OSQP, LinearAlgebra # temp
import controlutils
cu = controlutils
includet("MassSpringDamper.jl")

# create an instance
m = MassSpringDamperModel(0, # ma
	0, # ka
	500, # mo
	100) # umax
ny, nu = cu.dims(m)
opt = cu.OptOptions(false, false, 0.1, 1, :none, 1e-8, false)
fdes = 0.1 # KHz
N = round(Int, 1/(opt.fixedδt*fdes)) # 1 period with dt=0.1 in createInitialTraj
param0 = [20.0, # τ1
100.0, # ko
100.0, # bo
0.0] # τ2

# Generate a reference traj
trajt, traj0orig, trajt = createInitialTraj(m, opt, N, fdes)
σomax = cu.limits(m)[1]

# Make traj satisfy dyn constraint with these params?
# traj0 = cu.fixTrajWithDynConst(m, opt, traj0orig, param0)
traj0 = traj0orig

POPTS = cu.ParamOptOpts(
	τinds=[1,4], 
	R=(zeros(2,2), 0, 1.0*I), 
	plimsL = [0.1, 0.1, 0.1, 0.0],
	plimsU = [1000.0, 1000.0, 1000.0, 100.0],
	uinfnorm = false
)

σamax = 0.3 # [mm] constant? for robobee actuators
# σamax = 100 # [mm] constant? test EM

"""One-off ID or opt"""
function opt1(traj, param, mode, scaleTraj, bkratio=1.0, τ21ratiolim=2.0; testAffine=false, testAfter=false)
	# A polytope constraint for the params: simple bo >= ko =>  ko - bo <= 0. Second, τ2 <= 2*τ1 => -2*τ1 + τ2 <= 0
	Cp = Float64[0  bkratio  -1  0;
		-τ21ratiolim  0  0  1]
	dp = zeros(2)
	param1, paramObj, traj1, unactErr, paramConstraint, s = cu.optAffine(m, opt, traj, param, POPTS, mode, σamax; test=testAffine, scaleTraj=scaleTraj, Cp=Cp, dp=dp, print_level=1, max_iter=4000)
	if testAfter
		cu.affineTest(m, opt, traj1, param1, POPTS)
	end
	println("bkratio = ", bkratio, ", τ21ratiolim = ", τ21ratiolim, " => ", param1')
	return traj1, param1, paramObj, paramConstraint, s
end

"""Debug components in a traj. Assumes traj, param are feasible together here."""
function debugComponentsPlot(traj, param)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    τ1, ko, bo, τ2 = param
	# opt1(traj, param, 1, 1.0; testAffine=true)
	post, velt, acct = refTraj(m, fdes)
	tvec = trajt[1:N]

	# Get the components
	yo, Hτ, Hio, Hia, Hstiffo, Hstiffa, Hdamp = cu.paramAffine(m, opt, traj, param, POPTS; debugComponents=true)
	pt0, Tnew = cu.getpt(m, param)
	inertialo = zeros(N)
	inertiala = similar(inertialo)
	stiffo = similar(inertialo)
	stiffa = similar(inertialo)
	damp = similar(inertialo)

	for k=1:N
		inertialo[k] = (Hτ(Hio(yo(k), yo(k+1)), yo(k)) * pt0)[1]
		inertiala[k] = (Hτ(Hia(yo(k), yo(k+1)), yo(k)) * pt0)[1]
		stiffo[k] = (Hτ(Hstiffo(yo(k)), yo(k)) * pt0)[1]
		stiffa[k] = (Hτ(Hstiffa(yo(k)), yo(k)) * pt0)[1]
		damp[k] = (Hτ(Hdamp(yo(k)), yo(k)) * pt0)[1]
	end

	function plotComponents(ylbl)
		pl = plot(tvec, inertialo + inertiala, linewidth=2, label="i", ylabel=ylbl, legend=:outertopright)
		plot!(pl, tvec, stiffo, linewidth=2, label="so")
		# plot!(pl, trajt, stiffa, linewidth=2, label="sa")
		plot!(pl, tvec, damp, linewidth=2, label="d")
		tot = inertialo+inertiala+stiffo+stiffa+damp
		plot!(pl, tvec, tot, linewidth=2, linestyle=:dash, label="tot")
		# test what I think they should be
		y2, T, τfun, τifun = cu.transmission(m, traj, param; o2a=true)
    	plot!(pl, tvec, δt*T*m.mo*acct.(tvec), color=:black, linestyle=:dash, label="m*a")
    	plot!(pl, tvec, δt*T*ko*post.(tvec), color=:black, linestyle=:dash, label="k*x")
    	plot!(pl, tvec, δt*T*bo*velt.(tvec), color=:black, linestyle=:dash, label="b*dx")

		pl2 = plot(tvec, tot, linewidth=2, label="actn", legend=:outertopright)
		plot!(pl2, tvec, traj1[(N+1)*ny+1:end]*δt, linewidth=2, linestyle=:dash, label="act0")
		plot!(pl2, tvec, damp, linewidth=2, label="damp")
		return pl, pl2
	end

	# pl1 = plotTrajs(m, opt, trajt, [param], [traj])
	pls, plcomp = plotComponents("c")

	return pls, plcomp
end

function idealparams(param)
	# try to predict ideal params
	τ1, ko, bo, τ2 = param
	kidl = m.mo * (fdes * 2 * π)^2
	return [τ1, kidl, kidl, τ2] # Note the ko<=bo constraint: this is reflecting that, but must be set here separately
end

function plotNonlinBenefit()
    # First plot the param landscape
    pranges = [
        0:0.15:3.0, # τ21ratiolim
        0.1:0.04:0.8 # bkratios
    ]
    labels = [
        "nonlin ratio",
        "bkratio"
	]
	
	function maxu(τ21ratiolim, bkratio)
		traj1, param1, paramObj, paramConstraint, s = opt1(traj0, param0, 1, 1.0, bkratio, τ21ratiolim)
		return norm(traj1[(N+1)*ny+1:end], Inf)
	end

	function plotSlice(i1, i2)
		zgrid = [maxu(x,y) for y in pranges[i2], x in pranges[i1]] # reversed: see https://github.com/jheinen/GR.jl/blob/master/src/GR.jl
		# to get the improvement, divide each metric by the performance at τ2=0
		maxuatτ2_0 = zgrid[:,1]
		zgrid = zgrid ./ repeat(maxuatτ2_0, 1, length(pranges[i1]))

		pl = contourf(pranges[i1], pranges[i2], zgrid, fill=true, seriescolor=cgrad(:bluesreds), xlabel=labels[i1], ylabel=labels[i2])
        # just in case
        xlims!(pl, (pranges[i1][1], pranges[i1][end]))
        ylims!(pl, (pranges[i2][1], pranges[i2][end]))
        return pl
    end
    
    return (plotSlice(1, 2),)
end

# Test feasibility
function paramTest(p, s)
	xtest = [p; zeros((N+1)*ny); s]
	g = paramConstraint(xtest)
	# no unact constraint. have [gpolycon (1); gtransmission (1); ginfnorm (2*N)]
	println("Feas: should be nonpos: ", maximum(g), "; transmission: ", g[2])
	unew = cu.getTrajU(m, opt, traj1, p, POPTS)
	println("Obj: ", paramObj(xtest), " should be ", norm(unew, Inf)^2)
end

# One-off ID or opt ---------

# first optimization to get better params - closer to resonance
traj1, param1, paramObj, paramConstraint, s = opt1(traj0, param0, 1, 1.0, 0.2)
display(param1')
# param1 = idealparams(param1)

# debug components ---
pls = debugComponentsPlot(traj1, param1)
plot(pls..., size=(800,300))

# pls = plotNonlinBenefit() # SLOW
# plot(pls...)

# 
# Δτ1s = collect(0:0.1:15)
# unorms = unormΔτ1.(Δτ1s)

# # many sims (scale) --------------

# function paramsFor(σamax, scaleTraj)
# 	Tmin = 10.08799170499444*scaleTraj/σamax
# 	plimsL = [Tmin, 0.1, 0.1]
# 	plimsU = [1000.0, 1000.0, 1000.0]
# 	param1, _, traj1, _ = cu.optAffine(m, opt, traj0, param0, 1, R_WTS, 0.1, plimsL, plimsU, scaleTraj; Fext_pdep=true, test=false, testTrajReconstruction=false, print_level=1, max_iter=100)
# 	return [param1; norm(traj1[(N+1)*ny:end], Inf)]
# end

# scales = 0.1:0.2:2.0
# σamaxs = 0.1:1:20
# llabels = [
# 	"T",
# 	"k",
# 	"b"
# ]

# res = hcat(paramsFor.(0.3, scales)...)'
# np = length(param0)
# p1 = plot(scales, res[:,1:np], xlabel="traj scale", label=llabels, ylabel="design params", linewidth=2, legend=:topleft)
# p2 = plot(scales, res[:,np+1], xlabel="traj scale", ylabel="umin [mN]", linewidth=2, legend=false)
# res2 = hcat(paramsFor.(σamaxs, 1.0)...)'
# p3 = plot(σamaxs, res2[:,1:np], xlabel="max stroke [mm]", label=llabels, ylabel="design params", linewidth=2, legend=:topleft)
# p4 = plot(σamaxs, res2[:,np+1], xlabel="max stroke [mm]", ylabel="umin [mN]", linewidth=2, legend=false)
# plot(p1, p2, p3, p4)
