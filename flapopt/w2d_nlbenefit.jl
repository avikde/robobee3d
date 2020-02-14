using Dierckx, Plots

NLBENEFIT_FNAME = "nonlin.zip"
function nonlinBenefit(fname, ret, Tratios, minals, Qdts=[5e3, 1e4], Phis=[90,120], kos=range(30,60,length=4), wdens=range(0.009, 0.015, length=4); τ2eq=false, kwargs...)
	i = 0
	Ntotal = length(Tratios)*length(minals)*length(Qdts)*length(Phis)*length(kos)*length(wdens)
	function maxu(τ21ratiolim, minal, Qdt, phi, ko, wdens1)
		i += 1
		print(i,"/",Ntotal,": ")
		m.kbo[1] = ko
		rr = opt1(m, ret["traj"], ret["param"], 1, minal, τ21ratiolim; τ2eq=τ2eq, tol=5e-2, Qdt=Qdt, Φ=phi, wingdens1=wdens1, kwargs...)
		if rr["status"] >= -2
			return [rr["u∞"]; rr["al"]; rr["FD∞"]; rr["param"]]
		else
			println("BAD will skip")
			return [NaN; NaN; NaN; NaN]
		end
	end

	res = [[T;al;Qdt;phi;ko; maxu(T,al,Qdt,phi,ko,wdens1)] for T in Tratios, al in minals, Qdt in Qdts, phi in Phis, ko in kos, wdens1 in wdens]
	matwrite(fname, Dict("res" => res); compress=true)
	return res
end

function plotNonlinBenefit(fname, ypl; s=100, xpl=[0,3])
	results = matread(fname)["res"]
	
	Nhead = 5#length(size(results))

	function normalizeTbenefit(results)
		# Row 1 has Tratio=0 (res[1,:])
		for Tratio=size(results,1):-1:1
			for minal=1:size(results,2)
				resT0 = results[1,minal]
				@assert !isnan(resT0[Nhead+1]) "Linear transmission result was infeas"
				# normalize
				results[Tratio,minal][Nhead+1] /= resT0[Nhead+1]
				# results[Tratio,minal][Nhead+3] /= resT0[Nhead+3]
			end
		end
		return results
	end

	MODE = 1

	if MODE == 0
		results = results[:,:,1,1,1]
		results = normalizeTbenefit(results)
		ylabel = "FL [mg]"
	elseif MODE == 1
		results = results[:,3,1,2,:,3] # Tratio,ko
		ylabel = "ko ratio"
	elseif MODE == 2
		results = results[2,3,1,2,:,:] # Tratio,wingdens
		ylabel = "Izz"
	end
	# println(size(results))

	xyzi = zeros(length(results[1,1]),0)
	for res in results
		if !isnan(res[Nhead+1])
			xyzi = hcat(xyzi, res)
		end
	end
	# for plotting
	head = xyzi[1:Nhead,:]
	stats = xyzi[Nhead+1:Nhead+4,:]
	params = xyzi[Nhead+4:end,:]
	AL = stats[2,:]
	Tbenefit = stats[1,:]#clamp.(stats[1,:], 0.0, 1.0)#spline fit cause > 1 which doesn't make sense
	FLspec = AL./stats[3,:]
	ko = head[5,:]
	mw = params[3,:]
	Aw = params[6,:]
	cb2 = params[1,:]
	Lw = Aw ./ sqrt.(cb2)
	Izz = mw .* Lw.^2
	ko_I = head[5,:]./Izz
	T1 = params[2,:]
	koratio = ko./(ko + m.ka./T1.^2)
	Tractual = params[5,:]./params[2,:]

	# SET AXES HERE

	function contourFromUnstructured(xi, yi, zi; title="", rev=false)
		xpl = minimum(xi), maximum(xi)
		ypl = minimum(yi), maximum(yi)
		X = range(xpl..., length=50)
		Y = range(ypl..., length=50)
		# Spline from unstructured data https://github.com/kbarbary/Dierckx.jl
		# println("Total points = ", length(xi))
		spl = Spline2D(xi, yi, zi; s=s)
		function ff(x,y)
			zspl = spl(x,y)
			return zspl >= minimum(zi) && zspl <= maximum(zi) ? zspl : NaN
		end
		return contour(X, Y, ff, 
			titlefontsize=10, grid=false, lw=2, c=(rev ? :bluesreds_r : :bluesreds), 
			xlabel="T ratio", ylabel=ylabel, title=title,
			xlims=xpl, ylims=ypl)
	end

	if MODE == 0
		return [
			# scatter(xyzi[1,:], xyzi[4,:]),
			# scatter3d(xyzi[1,:], xyzi[4,:], xyzi[3,:]),
			contourFromUnstructured(Tractual, AL, Tbenefit; title="Nonlinear transmission benefit [ ]"),
			contourFromUnstructured(Tractual, AL, FLspec; title="FL sp.", rev=true),
			# contourFromUnstructured(Tractual, AL, params[2,:]; title="T1 [rad/mm]"),
			contourFromUnstructured(Tractual, AL, params[6,:]; title="Aw [mm^2]"),
			contourFromUnstructured(Tractual, AL, 1000.0 ./(N*params[7,:]); title="Freq [Hz]")
		]
	elseif MODE == 1
		return [
			contourFromUnstructured(Tractual, koratio, FLspec; title="FL sp.", rev=true),
			contourFromUnstructured(Tractual, koratio, ko_I; title="ko/I"),
			contourFromUnstructured(Tractual, koratio, params[6,:]; title="Aw [mm^2]"),
			contourFromUnstructured(Tractual, koratio, 1000.0 ./(N*params[7,:]); title="Freq [Hz]")
		]
	elseif MODE == 2
		return [
			scatter(Tractual, Izz),
			scatter3d(Tractual, Izz, FLspec)
			# contourFromUnstructured(Tractual, Izz, FLspec; title="FL sp.", rev=true)#,
			# contourFromUnstructured(Tractual, Izz, ko_I; title="ko/I"),
			# contourFromUnstructured(Tractual, Izz, params[6,:]; title="Aw [mm^2]"),
			# contourFromUnstructured(Tractual, Izz, 1000.0 ./(N*params[7,:]); title="Freq [Hz]")
		]
	end
end
