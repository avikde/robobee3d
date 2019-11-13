
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
using SparseArrays, OSQP, LinearAlgebra # temp
import controlutils
cu = controlutils
includet("MassSpringDamper.jl")

# create an instance
m = MassSpringDamperModel(6, # ma
	150, # ka
	500, # mo
	100) # umax
ny, nu = cu.dims(m)
opt = cu.OptOptions(false, 0.1, 1, :none, 1e-8, false)
fdes = 0.1 # KHz
N = round(Int, 10/fdes) # 1 period with dt=0.1 in createInitialTraj
param0 = [20.0, # τ1
100.0, # ko
10.0, # bo
20.0] # τ2

# Using the plant model for now, but eventually the ref traj should be generated by a template? TODO:
trajt, traj0orig, trajt = createInitialTraj(m, opt, N, fdes, [1e6, 1e5], param0, 227; showPlot=false)

# FIXME: replace this with a getOutputTrajectory(), and get σomax from that
const σomax = 9.970576473992493

# Make traj satisfy dyn constraint with these params?
traj0 = cu.fixTrajWithDynConst(m, opt, traj0orig, param0)

R_WTS = (zeros(2,2), 0, 1.0*I)#diagm(0=>[0.1,100]))
σamax = 0.3 # [mm] constant? for robobee actuators
# σamax = 100 # [mm] constant? test EM
Tmin = 10.08799170499444/σamax
plimsL = [0.1, 0.1, 0.1, 0.0]
plimsU = [1000.0, 1000.0, 1000.0, 100.0]

"""One-off ID or opt"""
function opt1(traj, param, mode, scaleTraj; testAffine=false, testAfter=false)
	optoptions = (R_WTS, 0.1, plimsL, plimsU, σamax)
	param1, paramObj, traj1, _ = cu.optAffine(m, opt, traj, param, mode, [1,4], optoptions..., scaleTraj, false; Fext_pdep=true, test=testAffine, testTrajReconstruction=false, print_level=1, max_iter=4000)
	if testAfter
		cu.optAffine(m, opt, traj1, param1, 2, optoptions...; Fext_pdep=true, test=true, testTrajReconstruction=false, print_level=1, max_iter=200) # TEST
	end
	return traj1, param1, paramObj
end

"""Debug components in a traj. Assumes traj, param are feasible together here."""
function debugComponentsPlot(traj, param)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)

	# Get the components
	yo, HMqTo, HMqTa, HCgJTo, HCgJTa = cu.paramAffine(m, opt, traj, param, R_WTS; Fext_pdep=true, debugComponents=true)
	pt0, Tnew = cu.getpt(m, param)
	inertialo = zeros(N)
	inertiala = similar(inertialo)
	stiffdampo = similar(inertialo)
	stiffdampa = similar(inertialo)

	for k=1:N
		# TODO: check *pt0
		inertialo[:,k] = (HMqTo(yo(k), yo(k+1)) - HMqTo(yo(k), yo(k))) * pt0
		inertiala[:,k] = (HMqTa(yo(k), yo(k+1)) - HMqTo(yo(k), yo(k))) * pt0
		stiffdampo[:,k] = (δt * HCgJTo(yo(k))) * pt0
		stiffdampa[:,k] = (δt * HCgJTa(yo(k))) * pt0
	end

	function plotComponents(ylbl)
		pl = plot(inertialo + inertiala, linewidth=2, label="i", ylabel=ylbl, legend=:outertopright)
		plot!(pl, stiffdampo, linewidth=2, label="g")
		plot!(pl, stiffdampa, linewidth=2, label="ga")
		tot = inertialo+inertiala+stiffdampo+stiffdampa
		plot!(pl, tot, linewidth=2, linestyle=:dash, label="tot")

		pl2 = plot(traj1[(N+1)*ny+1:end], linewidth=2, label="actf", legend=:outertopright)
		return pl, pl2
	end

	pl1 = plotTrajs(m, opt, trajt, [param], [traj])
	pls, plcomp = plotComponents("c")

	return pl1..., pls, plcomp
end

# One-off ID or opt ---------

traj1, param1, _ = opt1(traj0, param0, 1, 1.0)
display(param1')

# debug components ---

# pls = debugComponentsPlot(traj1, param1)
# plot(pls..., size=(800,600))
# gui()
# error("comp")

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
