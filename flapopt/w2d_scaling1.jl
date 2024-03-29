
using MAT
include("w2d_model.jl")
include("w2d_paramopt.jl")

const SCALING1_FNAME = "scaling1.zip"
const SCALING2_FNAME = "scaling2.zip"

function scaling1(m::Wing2DOFModel, opt, traj, param, τ21ratiolim, xs, minlifts, Qdts; kwargs...)
	np = length(param)
	i = 0
	Ntotal = length(xs)*length(minlifts)*length(Qdts)
	function scaling1single(x, minlift, Qdt)
		i += 1
		print(i, "/", Ntotal, ": ")
		r = opt1(m, traj, param, 1, minlift, τ21ratiolim; Φ=x, tol=5e-3, Qdt=Qdt, kwargs...)
		return [x; minlift; Qdt; r["param"]; r["u∞"]; r["al"]; r["δact"]; mean(abs.(r["mechPow"])); r["FD∞"]]
	end
	results = [scaling1single(x, minlift, Qdt) for minlift in minlifts, x in xs, Qdt in Qdts]
	resdict = Dict(
		"results" => results
	)
	# ^ returns a 2D array result arrays
	matwrite(SCALING1_FNAME, resdict; compress=true)

	return resdict
end

function scaling1disp(resarg; useFDasFact=true, scatterOnly=false, xpl=nothing, ypl=nothing, s=nothing, Fnom=75, mactline=7000)
	np = length(param0)
	resdict = typeof(resarg) == String ? matread(resarg) : resarg
	mactRobobee = Fnom*σamax
	Nax = 3 #2

	# Produced unstructured xi,yi,zi data
	xi = Float64[]
	ARi = Float64[]
	Qdti = Float64[]
	Awi = Float64[]
	FLi = Float64[]
	Phii = Float64[]
	mli = Float64[]
	Lwi = Float64[]
	macti = Float64[]  # times Robobee act
	powi = Float64[]
	freqi = Float64[]
	Ti = Float64[]
	Tri = Float64[]
	wi = Float64[]
	for res in resdict["results"]
		Phi = deg2rad(res[1])
		param = res[Nax+1:Nax+np]
		stats = res[Nax+np+1:end]
		Lw = param[6]/sqrt(param[1])

		# Append to unstructured data
		append!(Phii, Phi)
		append!(mli, res[2])
		append!(Qdti, res[3])
		append!(Lwi, Lw)
		append!(Awi, param[6])
		append!(ARi, Lw/sqrt(param[1]))
		append!(xi, Phi*Lw)
		append!(FLi, stats[2])
		append!(macti, stats[3] * (useFDasFact ? stats[5] : stats[1])/mactRobobee)
		append!(powi, stats[4])
		append!(freqi, 1000/(N*param[np]))
		append!(Ti, param[2])
		append!(wi, param[4])
		append!(Tri, param[5]/param[2])
	end

	Phis = unique(Phii)
	minlifts = unique(mli)
	Qdts = unique(Qdti)
	println("Phis = ", Phis, " minlifts = ", unique(mli), " Qdts = ", Qdts)
	# res2 = filter(x -> (deg2rad(x[1]) ≈ Phis[3] && x[2] ≈ minlifts[3]), resdict["results"])
	# display(hcat(res2...))

	# Output the plots
	pl1 = scatter(xlabel="Phi", ylabel="Lw", legend=false)
	pl2 = scatter(xlabel="x", ylabel="FL", legend=false)
	scatter!(pl1, Phii, Lwi)
	scatter!(pl2, xi, FLi)

	retpl = [pl1, pl2]

	if !scatterOnly
		# Plot range changes depending on opt results (see scatter)
		if isnothing(xpl)
			xpl = [minimum(xi), maximum(xi)]
			ypl = [minimum(yi), maximum(yi)]
		end
		X = range(xpl[1], xpl[2], length=50)
		Y = range(ypl[1], ypl[2], length=50)
		if isnothing(s)
			s = length(xi)
		end

		"Spline from unstructured data https://github.com/kbarbary/Dierckx.jl"
		function splineFromUnstructured(xi, yi, zi; Qdtsi=nothing, kwargs...)
			if !isnothing(Qdtsi)
				ii = findall(x -> x≈Qdts[Qdtsi], Qdti)
				return Spline2D(xi[ii], yi[ii], zi[ii]; s=s)
			else
				return Spline2D(xi, yi, zi; s=s)
			end
		end

		function contourFromUnstructured(xi, yi, zi; Qdtsi=nothing, title="")
			spl = splineFromUnstructured(xi, yi, zi; Qdtsi=Qdtsi)
			function ff(x,y)
				zspl = spl(x,y)
				return zspl >= minimum(zi) && zspl <= maximum(zi) ? zspl : NaN
			end
			return contour(X, Y, ff, 
				titlefontsize=10, grid=false, lw=2, c=:bluesreds, 
				xlabel="x [mm]", ylabel="FL [mg]", title=title,
				xlims=xpl, ylims=ypl)
		end

		function contourFromUnstructured2(xi, yi, zi; Qdtsi=nothing, title="")
			spl = splineFromUnstructured(xi, yi, zi; Qdtsi=Qdtsi)
			function ff(x,y)
				zspl = spl(x,y)
				return zspl >= minimum(zi) && zspl <= maximum(zi) ? zspl : NaN
			end
			X = range(1.5, 2.5, length=50)
			Y = range(20, 35, length=50)
			return contour(X, Y, ff, 
				titlefontsize=10, grid=false, lw=2, c=:bluesreds, 
				xlabel="mact [x]", ylabel="pow [mW]", title=title,
				xlims=[1.5, 2.5], ylims=[20, 35])
		end

		"Get an isoline along an mact https://github.com/avikde/robobee3d/pull/140"
		function isoline!(pl, xi, yi, zi, mact, Qdtsi; kwargs...)
			spl = splineFromUnstructured(xi, yi, zi; Qdtsi=Qdtsi)
			plot!(pl, X, spl.(X, mact./X); kwargs...)
		end
		function plisolines!(pl1, pl2, pl3, Qdtsi; kwargs...)
			isoline!(pl1, xi, FLi, FLi, mactline, Qdtsi; label="FL", kwargs...)
			isoline!(pl1, xi, FLi, 10*powi, mactline, Qdtsi; label="10pow", kwargs...)
			isoline!(pl2, xi, FLi, 100*Ti, mactline, Qdtsi; label="100T", kwargs...)
			isoline!(pl2, xi, FLi, freqi, mactline, Qdtsi; label="f", kwargs...)
			isoline!(pl3, xi, FLi, Awi, mactline, Qdtsi; label="Aw", kwargs...)
			# isoline!(pl3, xi, FLi, 100*wi, mactline, Qdtsi; label="10*whinge", kwargs...)
			# isoline!(pl3, xi, FLi, 50*ARi, mactline, Qdtsi; label="50AR")
		end
		function plcontours(Qdtsi)
			plmact = contourFromUnstructured(xi, FLi, macti; Qdtsi=Qdtsi, title="mact [x Robobee]")
			plot!(plmact, X, mactline./X, lw=2, color=:black, ls=:dash, label="")
			return [plmact, 
			# scatter3d(xi, FLi, macti, camera=(90,40)),
			contourFromUnstructured(xi, FLi, powi; Qdtsi=Qdtsi, title="Avg mech pow [mW]"),
			# contourFromUnstructured(xi, FLi, Tri; Qdtsi=4, title="T ratio 2"),
			# # contourFromUnstructured(xi, FLi, Tri; Qdtsi=2, title="T ratio 2"),
			# # contourFromUnstructured(xi, FLi, Tri; Qdtsi=3, title="T ratio 3"),
			# # contourFromUnstructured(xi, FLi, Tri; Qdtsi=4, title="T ratio 4"),
			# contourFromUnstructured(xi, FLi, rad2deg.(Phii); Qdtsi=4, title="Phi"),
			# contourFromUnstructured(xi, FLi, Awi; Qdtsi=Qdtsi, title="Aw [mm^2]"),
			# contourFromUnstructured(xi, FLi, ARi; Qdtsi=Qdtsi, title="ARi"),
			# scatter3d(xi, FLi, powi, camera=(10,40)),
			# contourFromUnstructured(xi, FLi, rad2deg.(Phii); title="Phi"), 
			# contourFromUnstructured(xi, FLi, mli; title="ml"), 
			# contourFromUnstructured(xi, FLi, freqi; Qdtsi=Qdtsi, title="freq [Hz]"), 
			contourFromUnstructured(xi, FLi, Ti; Qdtsi=Qdtsi, title="T1 [rad/mm]"), 
			contourFromUnstructured2(macti, powi, FLi./macti; Qdtsi=Qdtsi, title="FL/mact"), 
			contourFromUnstructured2(macti, powi, xi; Qdtsi=Qdtsi, title="x [mm]"), 
			contourFromUnstructured2(macti, powi, Awi; Qdtsi=Qdtsi, title="Aw [mm^2]"), 
			contourFromUnstructured2(macti, powi, Ti; Qdtsi=Qdtsi, title="T [rad/mm]"), 
			contourFromUnstructured2(macti, powi, Tri; Qdtsi=Qdtsi, title="T ratio [ ]")
			]
		end
			
		pliso = [plot(xlabel="x [mm]", legend=:outertopright), plot(xlabel="x [mm]", legend=:outertopright), plot(xlabel="x [mm]", legend=:outertopright)]
		plisolines!(pliso..., 2)
		plisolines!(pliso..., 1, ls=:dash)
		# plisolines!(pliso..., 2, ls=:dashdot)
		plisolines!(pliso..., 4, ls=:dot)

		# pl1 = plot(xs, [res[6]/res[1] for res in results], xlabel="Phi", ylabel="Lw", lw=2)

		append!(retpl, plcontours(2))
		# append!(retpl, pliso)
	end
	
	return retpl
end

function scaling2(m::Wing2DOFModel, opt, traj, param, phi, Qdts, minlift, τ21ratiolim; kwargs...)
	np = length(param)
	function scaling2single(#= phi,  =#Qdt)
		r = opt1(m, traj, param, 1, minlift, τ21ratiolim; Φ=phi, Qdt=Qdt, kwargs...)
		return [phi; Qdt; r["param"]; r["u∞"]; r["al"]; r["δact"]; mean(abs.(r["mechPow"])); r["FD∞"]]
	end
	results = [scaling2single(#= phi,  =#Qdt) for Qdt in Qdts]
	resdict = Dict(
		"results" => results
	)
	# ^ returns a 2D array result arrays
	matwrite(SCALING2_FNAME, resdict; compress=true)

	return resdict
end


"""Run many opts to get the best params for a desired min lift"""
function scaleParamsForlift(ret, minlifts, τ21ratiolim; kwargs...)
	traj, param = ret["traj"], ret["param"]
	function maxuForMinAvgLift(al)
		r = opt1(m, traj, param, 1, al, τ21ratiolim; kwargs...)
		# kΨ, bΨ = param2[4:5]
		uu = r["traj"][(N+1)*ny:end]
		return [r["param"]; norm(uu, Inf); norm(uu, 2)/N; r["al"]]
	end
	llabels = [
		"chord",
		"T1",
		"mwing",
		"hinge k",
		"hinge b",
		"T2",
		"Aw",
		"dt"
	]
	minliftsmg = minlifts

	res = hcat(maxuForMinAvgLift.(minlifts)...)'
	display(res)
	np = length(param0)
	actualliftsmg = res[:,np+4]
	p1 = plot(actualliftsmg, res[:,POPTS.τinds], xlabel="avg lift [mg]", label=llabels[POPTS.τinds], ylabel="T1,T2", linewidth=2, legend=:topleft)
	p2 = plot(actualliftsmg, [res[:,np+1]  res[:,np+3]], xlabel="avg lift [mg]", ylabel="umin [mN]", linewidth=2, legend=:topleft, label=["inf","2","al"])
	hline!(p2, [75], linestyle=:dash, color=:black, label="robobee act")
	p3 = plot(actualliftsmg, 1000 ./ (N*res[:,np]), xlabel="avg lift [mg]", ylabel="Cycle freq [Hz]", linewidth=2, legend=false)
	return p1, p2, p3
end

