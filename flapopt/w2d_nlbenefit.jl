using Dierckx, Plots

NLBENEFIT_FNAME = "nonlin.zip"
function nonlinBenefit(fname, ret, Tratios, minals; τ2eq=false, kwargs...)
	i = 0
	Ntotal = length(Tratios)*length(minals)
	function maxu(τ21ratiolim, minal,Qdt,phi,ko)
		i += 1
		print(i,"/",Ntotal,": ")
		m.kbo[1] = ko
		rr = opt1(m, ret["traj"], ret["param"], 1, minal, τ21ratiolim; τ2eq=τ2eq, tol=5e-2, Qdt=Qdt,Φ=phi, kwargs...)
		if rr["status"] >= -2
			return [rr["u∞"]; rr["al"]; rr["FD∞"]; rr["param"]]
		else
			println("BAD will skip")
			return [NaN; NaN; NaN; NaN]
		end
	end

	res = [[T;al;Qdt;phi;ko; maxu(T,al,Qdt,phi,ko)] for T in Tratios, al in minals, Qdt in [5e3, 1e4], phi in [90,120], ko in range(30,60,length=5)]
	matwrite(fname, Dict("res" => res); compress=true)
	return res
end

function plotNonlinBenefit(fname, ypl; s=100, xpl=[0,3])
	results = matread(fname)["res"]
	
	println(size(results))
	# Row 1 has Tratio=0 (res[1,:])
	for r=size(results,1):-1:1
		for c=1:size(results,2)
			resT0 = results[1,c]
			@assert !isnan(resT0[3]) "Linear transmission result was infeas"
			# normalize
			results[r,c][3] /= resT0[3]
			results[r,c][5] /= resT0[5]
		end
	end
	xyzi = zeros(length(results[1,1]),0)
	for res in results
		if !isnan(res[3])
			xyzi = hcat(xyzi, res)
		end
	end
	# lift to mg
	params = xyzi[6:end,:]

	X = range(xpl[1], xpl[2], length=50)
	Y = range(ypl[1], ypl[2], length=50)

	function contourFromUnstructured(xi, yi, zi; title="")
		# Spline from unstructured data https://github.com/kbarbary/Dierckx.jl
		# println("Total points = ", length(xi))
		spl = Spline2D(xi, yi, zi; s=s)
		function ff(x,y)
			zspl = spl(x,y)
			return zspl >= minimum(zi) && zspl <= maximum(zi) ? zspl : NaN
		end
		return contour(X, Y, ff, 
			titlefontsize=10, grid=false, lw=2, c=:bluesreds, 
			xlabel="T ratio", ylabel="FL [mg]", title=title,
			xlims=xpl, ylims=ypl)
	end

	Tractual = params[5,:]./params[2,:]
    
	return [
		# scatter(xyzi[1,:], xyzi[4,:]),
		# scatter3d(xyzi[1,:], xyzi[4,:], xyzi[3,:]),
		contourFromUnstructured(Tractual, xyzi[4,:], #= clamp.( =#xyzi[3,:]#= , 0.0, 1.0) =#; title="Nonlinear transmission benefit [ ]"),# opt errors cause > 1 which doesn't make sense
		contourFromUnstructured(Tractual, xyzi[4,:], params[2,:]; title="T1 [rad/mm]"),
		contourFromUnstructured(Tractual, xyzi[4,:], params[6,:]; title="Aw [mm^2]"),
		# contourFromUnstructured(xyzi[1,:], xyzi[4,:], Tractual; title="T ratio")
		contourFromUnstructured(Tractual, xyzi[4,:], 1000.0 ./(N*params[7,:]); title="Freq [Hz]")
	]
end
