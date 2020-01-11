
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
	2000, # mo
	100) # umax
ny, nu = cu.dims(m)
opt = cu.OptOptions(false, false, 0.1, 1, :none, 1e-8, false)
fdes = 0.1 # KHz
N = round(Int, 1/(opt.fixedδt*fdes)) # 1 period with dt=0.1 in createInitialTraj
param0 = [20.0, # τ1
100.0, # ko
100.0, # bo
0.0, # τ2
0.1] # dt

# Generate a reference traj
trajt, traj0orig, trajt = createInitialTraj(m, opt, N, fdes)
σomax = cu.limits(m)[1]

# Make traj satisfy dyn constraint with these params?
# traj0 = cu.fixTrajWithDynConst(m, opt, traj0orig, param0)
traj0 = traj0orig

POPTS = cu.ParamOptOpts(
	τinds=[1,4], 
	R=(zeros(2,2), 0, 1.0*I), 
	plimsL = [0.1, 0.1, 0.1, 0.0, 0.01],
	plimsU = [1000.0, 1000.0, 1000.0, 100.0, 0.5],
	uinfnorm = false
)

σamax = 0.3 # [mm] constant? for robobee actuators
# σamax = 100 # [mm] constant? test EM

"""One-off ID or opt"""
function opt1(traj, param, mode, scaleTraj, bkratio=1.0, τ21ratiolim=2.0; testAffine=false, testAfter=false)
	# A polytope constraint for the params: simple bo >= bkratio*ko =>  bkratio*ko - bo <= 0. Second, τ2 <= 2*τ1 => -2*τ1 + τ2 <= 0
	Cp = Float64[0  bkratio  -1  0  0;
		-τ21ratiolim  0  0  1  0]
	dp = zeros(2)
	ret = cu.optAffine(m, opt, traj, param, POPTS, mode, σamax; test=testAffine, scaleTraj=scaleTraj, Cp=Cp, dp=dp, print_level=1, max_iter=4000)

	uu = ret["traj"][(N+1)*ny:end]
	ret["u∞"] = norm(uu, Inf)
	if testAfter
		cu.affineTest(m, opt, ret["traj"], ret["param"], POPTS)
	end

	println("st ", scaleTraj, ": ", ret["status"], ", ", round.(ret["param"]', digits=3), ", fHz=", round(1000/(N*ret["param"][end]), digits=1), ", u∞=", round(ret["u∞"], digits=1), ", J=", round(ret["eval_f"](ret["x"]), digits=1))
	return ret
end

"""Debug components in a traj. Assumes traj, param are feasible together here."""
function debugComponentsPlot(traj, param; refScale=1.0)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    τ1, ko, bo, τ2, dt = param
	tvec = dt*collect(0:(N-1)) # was trajt[1:N]
	freq = 1/(N*dt)
	post, velt, acct = refTraj(m, freq)

	# Get the components
	yo, Hio, Hia, Hstiffo, Hstiffa, Hdamp = cu.paramAffine(m, opt, traj, param, POPTS; debugComponents=true) # NOTE: scaleTraj=1 here
	pt0, Tnew = cu.getpt(m, param)
	inertialo = zeros(N)
	inertiala = similar(inertialo)
	stiffo = similar(inertialo)
	stiffa = similar(inertialo)
	damp = similar(inertialo)
	# to fill in for (non) inertial half of H https://github.com/avikde/robobee3d/pull/102
	H0 = zeros(1,length(pt0)÷2)

	# println("Mass, refl=", m.mo, ",", m.ma/τ1^2)
	# TODO: update this
	meq = m.mo + m.ma/τ1^2
	keq = ko + m.ka/τ1^2

	for k=1:N
		σo = yo(k)[1]
		inertialo[k] = [cu.Hτ(Hio(yo(k), yo(k+1)), σo)*dt  H0] ⋅ pt0
		inertiala[k] = [cu.Hτ(Hia(yo(k), yo(k+1)), σo)*dt  H0] ⋅ pt0
		stiffo[k] = [H0  cu.Hτ(Hstiffo(yo(k)), σo)] ⋅ pt0
		stiffa[k] = [H0  cu.Hτ(Hstiffa(yo(k)), σo)] ⋅ pt0
		damp[k] = [H0  cu.Hτ(Hdamp(yo(k)), σo)] ⋅ pt0
	end

	function plotComponents(ylbl)
		pl = plot(tvec, (inertialo + inertiala)/dt, linewidth=2, label="i", ylabel=ylbl, legend=:outertopright)
		plot!(pl, tvec, (stiffo+stiffa)/dt, linewidth=2, label="s")
		# plot!(pl, trajt, stiffa, linewidth=2, label="sa")
		plot!(pl, tvec, damp/dt, linewidth=2, label="d")
		tot = inertialo+inertiala+stiffo+stiffa+damp
		plot!(pl, tvec, tot/dt, linewidth=2, linestyle=:dash, label="tot")
		# test what I think they should be
		y2, T, τfun, τifun = cu.transmission(m, traj, param; o2a=true)
    	plot!(pl, tvec, refScale*T*meq*acct.(tvec), color=:blue, linestyle=:dash, lw=2, label="m*a")
    	plot!(pl, tvec, refScale*T*keq*post.(tvec), color=:red, linestyle=:dash, lw=2, label="k*x")
    	plot!(pl, tvec, refScale*T*bo*velt.(tvec), color=:green, linestyle=:dash, lw=2, label="b*dx")

		pl2 = plot(tvec, tot/dt, linewidth=2, label="act(tot)", legend=:outertopright)
		plot!(pl2, tvec, ret1["traj"][(N+1)*ny+1:end], linewidth=2, linestyle=:dash, label="act")
		plot!(pl2, tvec, damp/dt, linewidth=2, label="damp")
		return pl, pl2
	end

	pl1 = plotTrajs(m, opt, trajt, [param], [traj]; fdes=freq, refScale=refScale)
	pls, plcomp = plotComponents("c")

	return pl1..., pls, plcomp
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

function actForceForHeightsPlot(bkratio)
	function costFor(height, dt)
		STROKE_AMP = 10 # [mm]
		G_CONST = 9.81e-3 # [mm/ms^2]
		loVel = sqrt(2*G_CONST*height)
		# loVel = STROKE_AMP*omega = STROKE_AMP*scale*2*pi/(N*dt)
		scaleTraj = loVel*dt*N/(STROKE_AMP*2*pi)
		# Tmin = 10.08799170499444*scaleTraj/σamax
		POPTS.plimsU[end] = dt
		ret1 = opt1(traj0, param0, 1, scaleTraj, bkratio)
		return ret1["u∞"]#[ret1["u∞"]; ret1["param"]]
	end

	heights = 20:2:40
	dtmaxs = 0.2:0.2:1.0

	pl = contour(heights, dtmaxs, costFor, fill=true, seriescolor=cgrad(:bluesreds), xlabel="Height [mm]", ylabel="dt")
	# # TODO: figure out how to draw two sets of contours in Julia
	# f2(a,b)=a/b
	# pl2 = contour(heights, dtmaxs, f2, seriescolor=cgrad(:bluesreds), fill=true, xlabel="Amplitude", ylabel="dt")
	plot(pl)
end

# One-off ID or opt ---------

POPTS.plimsU[end] = 0.4
trajScale = 0.56
ret1 = opt1(traj0, param0, 1, trajScale, 0.1)
display(ret1["param"]')

# debug components ---
pls = debugComponentsPlot(ret1["traj"], ret1["param"]; refScale=trajScale)
plot(pls..., size=(800,400))

# # ---------
# actForceForHeightsPlot(0.1)

# pls = plotNonlinBenefit() # SLOW
# plot(pls...)


# Δτ1s = collect(0:0.1:15)
# unorms = unormΔτ1.(Δτ1s)

# # many sims (scale) --------------

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
