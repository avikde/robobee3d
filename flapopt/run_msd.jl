
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
function opt1(traj, param, mode, scaleTraj; testAffine=false, testAfter=false)
	optoptions = (R_WTS, 0.1, plimsL, plimsU, σamax)
	# A polytope constraint for the params: simple bo >= ko =>  ko - bo <= 0 => 
	Cp = Float64[0  1  -1  0]
	dp = [0.0]
	param1, paramObj, traj1, _ = cu.optAffine(m, opt, traj, param, mode, [1,4], optoptions..., scaleTraj, false, Cp, dp; Fext_pdep=true, test=testAffine, testTrajReconstruction=false, print_level=1, max_iter=4000)
	if testAfter
		cu.optAffine(m, opt, traj1, param1, 2, optoptions...; Fext_pdep=true, test=true, testTrajReconstruction=false, print_level=1, max_iter=200) # TEST
	end
	return traj1, param1, paramObj
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
		# FIXME: this is OK for τ2=0 but need to use transmission()
		inertialo[k] = (Hτ(Hio(yo(k), yo(k+1)), yo(k)) * pt0/τ1)[1]
		inertiala[k] = (Hτ(Hia(yo(k), yo(k+1)), yo(k)) * pt0/τ1)[1]
		stiffo[k] = (Hτ(Hstiffo(yo(k)), yo(k)) * pt0/τ1)[1]
		stiffa[k] = (Hτ(Hstiffa(yo(k)), yo(k)) * pt0/τ1)[1]
		damp[k] = (Hτ(Hdamp(yo(k)), yo(k)) * pt0/τ1)[1]
	end

	function plotComponents(ylbl)
		pl = plot(tvec, inertialo + inertiala, linewidth=2, label="i", ylabel=ylbl, legend=:outertopright)
		plot!(pl, tvec, stiffo, linewidth=2, label="so")
		# plot!(pl, trajt, stiffa, linewidth=2, label="sa")
		plot!(pl, tvec, damp, linewidth=2, label="d")
		tot = inertialo+inertiala+stiffo+stiffa+damp
		plot!(pl, tvec, tot, linewidth=2, linestyle=:dash, label="tot")
		# test what I think they should be
    	plot!(pl, tvec, δt*m.mo*acct.(tvec), color=:black, linestyle=:dash, label="m*a")
    	plot!(pl, tvec, δt*ko*post.(tvec), color=:black, linestyle=:dash, label="k*x")
    	plot!(pl, tvec, δt*bo*velt.(tvec), color=:black, linestyle=:dash, label="b*dx")

		pl2 = plot(tvec, traj1[(N+1)*ny+1:end]*δt/τ1, linewidth=2, label="actf", legend=:outertopright)
		plot!(pl2, tvec, tot, linewidth=2, linestyle=:dash, label="tot")
		plot!(pl2, tvec, damp, linewidth=2, label="d")
		return pl, pl2
	end

	pl1 = plotTrajs(m, opt, trajt, [param], [traj])
	pls, plcomp = plotComponents("c")

	return pl1..., pls, plcomp
end

# One-off ID or opt ---------

# first optimization to get better params - closer to resonance
traj1, param1, _ = opt1(traj0, param0, 1, 1.0)
display(param1')
# try to predict ideal params
τ1, ko, bo, τ2 = param1
paramidl = [τ1, m.mo * (fdes * 2 * π)^2, bo, τ2]
display(paramidl')

# debug components ---

pls = debugComponentsPlot(traj1, paramidl)
plot(pls..., size=(800,600))
gui()
error("comp")

# TEST manual params
Hk, yo, umeas, B, N = cu.paramAffine(m, opt, traj1, param1, R_WTS, 1.0; Fext_pdep=true)
Δy0 = zeros((N+1)*ny)
function testp(pnew)
	trajnew = cu.reconstructTrajFromΔy(m, opt, traj1, yo, Hk, B, Δy0, pnew, false)
	uvec = trajnew[(N+1)*ny+1:end]
	println(pnew[end], " RMS u = ", norm(uvec)/N, ", const = ", cu.gtransmission(m, pnew, σomax) - σamax)
	return trajnew, pnew, norm(uvec)/N
end
pp = copy(param1)
ppfeas(Δτ1) = pp + [-Δτ1, 0, 0, 3/σamax^2*Δτ1]
testΔτ1(Δτ1) = testp(ppfeas(Δτ1))[1:2]
traj2, p2 = testΔτ1(0)
traj3, p3 = testΔτ1(1)

unormΔτ1(Δτ1) = testp(ppfeas(Δτ1))[3]
Δτ1s = collect(0:0.1:15)
unorms = unormΔτ1.(Δτ1s)
pl2 = plot(Δτ1s, unorms, xlabel="Δτ1", ylabel="unorm")

pl1 = plotTrajs(m, opt, trajt, [param1, p2, p3], [traj1, traj2, traj3]; ulim=1e4)
plot(pl1..., pl2)


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
