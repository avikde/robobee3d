
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
opt = cu.OptOptions(false, 0.1, 1, :none, 1e-8, false)
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

R_WTS = (zeros(2,2), 0, 1.0*I)#diagm(0=>[0.1,100]))
σamax = 0.3 # [mm] constant? for robobee actuators
# σamax = 100 # [mm] constant? test EM
plimsL = [0.1, 0.1, 0.1, 0.0]
plimsU = [1000.0, 1000.0, 1000.0, 100.0]

"""One-off ID or opt"""
function opt1(traj, param, mode, scaleTraj, bkratio=1.0; testAffine=false, testAfter=false)
	println("bkratio = ", bkratio)
	optoptions = (R_WTS, 0.1, plimsL, plimsU, σamax)
	# A polytope constraint for the params: simple bo >= ko =>  ko - bo <= 0 => 
	Cp = Float64[0  bkratio  -1  0]
	dp = [0.0]
	param1, paramObj, traj1, unactErr, paramConstraint = cu.optAffine(m, opt, traj, param, mode, [1,4], optoptions..., scaleTraj, false, Cp, dp; Fext_pdep=true, test=testAffine, testTrajReconstruction=false, print_level=1, max_iter=4000)
	if testAfter
		cu.optAffine(m, opt, traj1, param1, 2, optoptions...; Fext_pdep=true, test=true, testTrajReconstruction=false, print_level=1, max_iter=200) # TEST
	end
	return traj1, param1, paramObj, paramConstraint
end

"""Debug components in a traj. Assumes traj, param are feasible together here."""
function debugComponentsPlot(traj, param)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    τ1, ko, bo, τ2 = param
	# opt1(traj, param, 1, 1.0; testAffine=true)
	post, velt, acct = refTraj(m, fdes)
	tvec = trajt[1:N]

	# Get the components
	yo, Hτ, Hio, Hia, Hstiffo, Hstiffa, Hdamp = cu.paramAffine(m, opt, traj, param, R_WTS; Fext_pdep=true, debugComponents=true)
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

function unormΔτ1(Δτ1, bkratio)
	traj1, param1, paramObj, xConstraint = opt1(traj0, param0, 1, 1.0, bkratio)
	pnew = copy(param1) + [-Δτ1, 0, 0, 3/σamax^2*Δτ1]
	Hk, yo, umeas, B, N = cu.paramAffine(m, opt, traj1, pnew, R_WTS, 1.0; Fext_pdep=true)	
	Δy0 = zeros((N+1)*ny)
	trajnew = cu.reconstructTrajFromΔy(m, opt, traj1, yo, Hk, B, Δy0, pnew, false)	
	return norm(trajnew[(N+1)*ny+1:end], Inf)
end

function plotNonlinBenefit()
    # First plot the param landscape
    pranges = [
        0:0.3:6.0, # Δτ1s
        0.1:0.05:1.0 # bkratios
    ]
    labels = [
        "Delta t1",
        "bkratio"
	]
	
	function unormimprovement(Δτ1, bkratio)
		r = unormΔτ1(Δτ1, bkratio)/unormΔτ1(0, bkratio) # FIXME: this runs the 0 one twice...
		return r > 1.1 ? NaN : r
	end

    function plotSlice(i1, i2)
		pl = contour(pranges[i1], pranges[i2], unormimprovement, fill=true, seriescolor=cgrad(:bluesreds), xlabel=labels[i1], ylabel=labels[i2])
		# f1(Δτ1) = unormΔτ1(Δτ1, 0.2)
		# yy = f1.(pranges[i1])
		# println("hi", yy)
		# pl = plot(pranges[i1], yy)
        # just in case
        xlims!(pl, (pranges[i1][1], pranges[i1][end]))
        ylims!(pl, (pranges[i2][1], pranges[i2][end]))
        return pl
    end
    
    return (plotSlice(1, 2),)
end

# One-off ID or opt ---------

# first optimization to get better params - closer to resonance
traj1, param1, paramObj, xConstraint = opt1(traj0, param0, 1, 1.0, 0.2)
# Test sufficient approx for feasible transmission params
pp = copy(param1)
ppfeas(Δτ1) = pp + [-Δτ1, 0, 0, 3/σamax^2*Δτ1]
# Test feasibility
gparam(p) = xConstraint([p; zeros((N+1)*ny)])
display(param1')
# param1 = idealparams(param1)

# debug components ---
pls = debugComponentsPlot(traj1, ppfeas(4))
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
