using Dierckx, Plots

NLBENEFIT_FNAME = "nonlin.zip"
function nonlinBenefit(fname, ret, Tratios, minals, Qdts=[5e3, 1e4], Phis=[90,120], kos=range(30,60,length=5); τ2eq=false, kwargs...)
	i = 0
	Ntotal = length(Tratios)*length(minals)*length(Qdts)*length(Phis)*length(kos)
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

	res = [[T;al;Qdt;phi;ko; maxu(T,al,Qdt,phi,ko)] for T in Tratios, al in minals, Qdt in Qdts, phi in Phis, ko in kos]
	matwrite(fname, Dict("res" => res); compress=true)
	return res
end

function plotNonlinBenefit(fname, ypl; s=100, xpl=[0,3])
	results = matread(fname)["res"]
	
	Nhead = length(size(results))
	results = results[:,:,1,1,1] # pick Qdt, phi
	println(size(results))
	# Row 1 has Tratio=0 (res[1,:])
	for Tratio=size(results,1):-1:1
		for minal=1:size(results,2)
			resT0 = results[1,minal]
			@assert !isnan(resT0[Nhead+1]) "Linear transmission result was infeas"
			# normalize
			results[Tratio,minal][Nhead+1] /= resT0[Nhead+1]
			results[Tratio,minal][Nhead+3] /= resT0[Nhead+3]
		end
	end
	xyzi = zeros(length(results[1,1]),0)
	for res in results
		if !isnan(res[Nhead+1])
			xyzi = hcat(xyzi, res)
		end
	end
	# for plotting
	stats = xyzi[Nhead+1:Nhead+4,:]
	params = xyzi[Nhead+4:end,:]
	AL = stats[2,:]
	Tbenefit = stats[1,:]#clamp.(stats[1,:], 0.0, 1.0)

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
		contourFromUnstructured(Tractual, AL, Tbenefit; title="Nonlinear transmission benefit [ ]"),# opt errors cause > 1 which doesn't make sense
		# contourFromUnstructured(Tractual, AL, params[2,:]; title="T1 [rad/mm]"),
		contourFromUnstructured(Tractual, AL, params[6,:]; title="Aw [mm^2]"),
		# contourFromUnstructured(xyzi[1,:], xyzi[4,:], Tractual; title="T ratio")
		contourFromUnstructured(Tractual, AL, 1000.0 ./(N*params[7,:]); title="Freq [Hz]")
	]
end
