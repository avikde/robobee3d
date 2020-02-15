using Dierckx, Plots

NLBENEFIT_FNAME = "nonlin.zip"
function nonlinBenefit(fname, ret, Tratios, minals, Qdts=[5e3, 1e4], Phis=[90,120], kos=range(25,65,length=6), wdens=range(0.009, 0.013, length=3); τ2eq=false, kwargs...)
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
	resultsOrig = matread(fname)["res"]
	
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

	function contourFromUnstructured(xi, yi, zi, ylabel; colorbar_title="", title="", rev=false, ypl=nothing)
		xpl = minimum(xi), maximum(xi)
		if isnothing(ypl)
			ypl = minimum(yi), maximum(yi)
		end
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
			xlabel="T ratio", ylabel=ylabel, colorbar_title=colorbar_title, title=title,
			xlims=xpl, ylims=ypl)
	end

	function createPlots(MODE, ii; title="", ypl=nothing)
		# T;al;Qdt;phi;ko
		if MODE == 0
			results = resultsOrig[:,:,1,1,1,1]
			results = normalizeTbenefit(results)
			ylabel = "FL [mg]"
		elseif MODE == 1
			results = resultsOrig[:,ii[1],ii[2],ii[3],:,ii[4]] # Tratio,ko
			ylabel = "ko ratio"
		elseif MODE == 2
			results = resultsOrig[ii[1],ii[2],ii[3],ii[4],:,:] # Tratio,wingdens
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
		# FLspec = AL./stats[3,:]
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

		if MODE == 0
			return [
				# scatter(xyzi[1,:], xyzi[4,:]),
				# scatter3d(xyzi[1,:], xyzi[4,:], xyzi[3,:]),
				contourFromUnstructured(Tractual, AL, Tbenefit, ylabel; title="Nonlinear transmission benefit [ ]"),
				contourFromUnstructured(Tractual, AL, FLspec, ylabel; title="FL sp.", rev=true),
				# contourFromUnstructured(Tractual, AL, params[2,:]; title="T1 [rad/mm]"),
				contourFromUnstructured(Tractual, AL, params[6,:], ylabel; title="Aw [mm^2]"),
				contourFromUnstructured(Tractual, AL, 1000.0 ./(N*params[7,:]), ylabel; title="Freq [Hz]")
			]
		elseif MODE == 1
			return [
				# scatter3d(Tractual, koratio, FLspec),
				contourFromUnstructured(Tractual, koratio, FLspec, ylabel; colorbar_title="FL sp.", title=title, rev=true, ypl=ypl),
				# contourFromUnstructured(Tractual, koratio, ko_I, ylabel; title="ko/I"),
				# contourFromUnstructured(Tractual, koratio, params[6,:], ylabel; title="Aw [mm^2]"),
				# contourFromUnstructured(Tractual, koratio, 1000.0 ./(N*params[7,:]), ylabel; title="Freq [Hz]")
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
	
	return [
		# al;Qdt;phi;wingdens
		plot(createPlots(1, [3, 1, 2, 2]; title="Low I, 120deg")..., createPlots(1, [3, 1, 2, 3]; title="High I, 120deg")..., layout=(2,1)),
		plot(createPlots(1, [2, 1, 1, 2]; title="Low I, 90deg", ypl=(0.5,0.625))..., createPlots(1, [2, 1, 1, 3]; title="High I, 90deg", ypl=(0.5,0.625))..., layout=(2,1)),
		# plot(createPlots(1, [4, 1, 2, 2])..., createPlots(1, [4, 1, 2, 3])..., title="Wing density (low, high)"),
	]
	# createPlots(1, [3, 1, 2, 3])
	# createPlots(2, [2,3,1,2])
end

function nlNormalizeByOutput(m, opt, param)
	# TODO:
	param[2] *= 0.95
	param[5] = 2*param[2]
	return param
end

"""Generate plot like in [Jafferis (2016)] Fig. 4"""
function openLoopPlot(m, opt, param0, Vmin, Vmax; save=false)
	function getResp(f, uamp, nlt=false)
		param = copy(param0)
		if nlt
			param = nlNormalizeByOutput(m, opt, param)
		end
		ts = createInitialTraj(m, opt, 0, f, [1e3, 1e2], param, 0; uampl=uamp, trajstats=true, thcoeff=0.1)
		println("act disp=",ts[end])
		return ts[1:2]
	end
	fs = 0.03:0.005:0.25
	mN_PER_V = 75/160

	p1 = plot(ylabel="Norm. stroke ampl [deg/V]", ylims=(0.3,0.8))
	p2 = plot(xlabel="Freq [kHz]", ylabel="Hinge ampl [deg]", legend=false, ylims=(0,100))

	function plotForTrans(nlt)
		nltstr = nlt ? "N" : "L"
		for Vamp=range(Vmin, Vmax,length=2)
			println("Openloop @ ", Vamp, "V ", nltstr)
			uamp = Vamp*mN_PER_V
			amps = hcat(getResp.(fs, uamp, nlt)...)
			amps *= 180/pi # to degrees
			amps[1,:] /= (Vamp) # normalize
			amps[2,:] /= 2.0 # hinge ampl one direction
			# println(amps)
			plot!(p1, fs, amps[1,:], lw=2, label=string(nltstr, Vamp,"V"), ls=nlt ? :solid : :dash)
			plot!(p2, fs, amps[2,:], lw=2, label=string(nltstr, Vamp,"V"), ls=nlt ? :solid : :dash)
		end
	end

	plotForTrans(false)
	plotForTrans(true)

	plot(p1, #= p2, layout=(2,1),  =# size=(400, 300), dpi=200)
	println("dens=", param0[3]/param0[6], ", koratio=", m.kbo[1]/(m.kbo[1] + m.ka/param0[2]^2))
	if save
		savefig("olplot.png")
	end
	gui()
	error("Open loop plot")
end
