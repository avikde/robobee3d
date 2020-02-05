
const SCALING1_FNAME = "scaling1.zip"
function scaling1(m::Wing2DOFModel, opt, traj, param, xs, minlifts, τ21ratiolim; kwargs...)
	np = length(param)
	function scaling1single(x, minlift)
		r = opt1(m, traj, param, 1, minlift, τ21ratiolim; Φ=x, kwargs...)
		return [x; minlift; r["param"]; r["u∞"]; r["al"]; r["δact"]; mean(abs.(r["mechPow"])); r["FD∞"]; norm(r["unactErr"], Inf)]
	end
	results = [scaling1single(x, minlift) for minlift in minlifts, x in xs] # reversed
	resdict = Dict(
		"Phis" => xs, "minlifts" => minlifts, "results" => results
	)
	# ^ returns a 2D array result arrays
	matwrite(SCALING1_FNAME, resdict; compress=true)

	return resdict
end

function scaling1disp(resarg; useFDasFact=true, scatterOnly=false, xpl=nothing, ypl=nothing, s=nothing, Fnom=75, mactline=7000)
	np = length(param0)
	resdict = typeof(resarg) == String ? matread(resarg) : resarg
	mactRobobee = Fnom*σamax

	# Produced unstructured xi,yi,zi data
	xi = Float64[]
	ARi = Float64[]
	Awi = Float64[]
	FLi = Float64[]
	Phii = Float64[]
	mli = Float64[]
	Lwi = Float64[]
	macti = Float64[]  # times Robobee act
	powi = Float64[]
	freqi = Float64[]
	Ti = Float64[]
	for res in resdict["results"]
		Phi = deg2rad(res[1])
		param = res[2+1:2+np]
		stats = res[2+np+1:end]
		Lw = param[6]/sqrt(param[1])

		# Append to unstructured data
		append!(Phii, Phi)
		append!(mli, res[2])
		append!(Lwi, Lw)
		append!(Awi, param[6])
		append!(ARi, Lw/sqrt(param[1]))
		append!(xi, Phi*Lw)
		append!(FLi, stats[2])
		append!(macti, stats[3] * (useFDasFact ? stats[5] : stats[1])/mactRobobee)
		append!(powi, stats[4])
		append!(freqi, 1000/(N*param[np]))
		append!(Ti, param[2])
	end

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

		function contourFromUnstructured(xi, yi, zi; title="")
			# Spline from unstructured data https://github.com/kbarbary/Dierckx.jl
			# println("Total points = ", length(xi))
			spl = Spline2D(xi, yi, zi; s=s)
			ff(x,y) = spl(x,y)
			return contour(X, Y, ff, 
				titlefontsize=10, grid=false, lw=2, c=:bluesreds, 
				xlabel="x [mm]", ylabel="FL [mg]", title=title,
				xlims=xpl, ylims=ypl)
		end

		# pl1 = plot(xs, [res[6]/res[1] for res in results], xlabel="Phi", ylabel="Lw", lw=2)

		plmact = contourFromUnstructured(xi, FLi, macti; title="mact [x Robobee]")
		plot!(plmact, X, mactline./X, lw=2, color=:black, ls=:dash, label="")

		append!(retpl, [plmact, 
			# scatter3d(xi, FLi, macti, camera=(90,40)),
			contourFromUnstructured(xi, FLi, powi; title="Avg mech pow [mW]"),
			contourFromUnstructured(xi, FLi, Awi; title="Aw [mm^2]"),
			contourFromUnstructured(xi, FLi, ARi; title="ARi"),
			# scatter3d(xi, FLi, powi, camera=(10,40)),
			# contourFromUnstructured(xi, FLi, rad2deg.(Phii); title="Phi"), 
			# contourFromUnstructured(xi, FLi, mli; title="ml"), 
			contourFromUnstructured(xi, FLi, freqi; title="freq [Hz]"), 
			contourFromUnstructured(xi, FLi, Ti; title="T1 [rad/mm]")])
	end
	
	return retpl
end


"""Run many opts to get the best params for a desired min lift"""
function scaleParamsForlift(ret, minlifts, τ21ratiolim; kwargs...)
	traj, param = ret["traj"], ret["param"]
	function maxuForMinAvgLift(al)
		r = opt1(m, traj, param, 1, al, τ21ratiolim; kwargs...)
		# kΨ, bΨ = param2[4:5]
		uu = r["traj"][(N+1)*ny:end]
		return [r["param"]; norm(uu, Inf); norm(r["unactErr"], Inf); norm(uu, 2)/N; r["al"]]
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
	p4 = plot(actualliftsmg, res[:,np+2], xlabel="avg lift [mg]", ylabel="unact err", linewidth=2, legend=false)

	return p1, p2, p3, p4
end
