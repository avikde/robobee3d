using Dierckx, Plots, MAT

NLBENEFIT_FNAME = "nonlin.zip"
function nonlinBenefit(fname, ret, Tratios, minals, Qdts=[5e3, 1e4], Phis=[90,120], kos=range(20,70,length=8), wdens=range(0.01, 0.013, length=2); τ2eq=false, kwargs...)
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

"Helper for plotNonlinBenefit"
function prepareDataNonlinBenefit(fname, MODE, ii)
	Nhead = 5#length(size(results))
	resultsOrig = matread(fname)["res"]

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
	
	# Create unstructured data array
	xyzi = zeros(length(first(resultsOrig)),0)
	function appendunstruc!(v, res)
		if !isnan(res[Nhead+1])
			v = hcat(v, res)
		end
		return v
	end

	# T;al;Qdt;phi;ko
	if MODE == 0
		results = resultsOrig[:,:,1,1,1,1]
		results = normalizeTbenefit(results)
		for res in results
			xyzi = appendunstruc!(xyzi, res)
		end
		ylabel = "FL [mg]"
	elseif MODE == 1
		for res in resultsOrig[:,ii[1],ii[2],ii[3],:,ii[4]] # Tratio,ko
			xyzi = appendunstruc!(xyzi, res)
		end
		ylabel = "ko ratio"
	elseif MODE == 2
		for res in resultsOrig[ii[1],ii[2],ii[3],ii[4],:,:] # Tratio,wingdens
			xyzi = appendunstruc!(xyzi, res)
		end
		ylabel = "Izz"
	end
	println("Data:", size(resultsOrig), "->", size(xyzi))
	return xyzi, Nhead, ylabel
end

"Plot benefits of nonlin transmission from big sim data"
function plotNonlinBenefit(fname, ypl; s=100, xpl=[0,3])
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

	function scatterAndShadow(xi, yi, zi, camera, ylabel; markeralpha=0.1)
		Z = ones(size(zi))
		p1 = scatter3d(camera=camera, legend=false, xlabel="T ratio", ylabel=ylabel)
		# shadows
		scatter3d!(p1, xi, yi, minimum(zi)*Z, color=:gray, markeralpha=markeralpha, markerstrokewidth=0)
		scatter3d!(p1, minimum(xi)*Z, yi, zi, color=:gray, markeralpha=markeralpha, markerstrokewidth=0)
		scatter3d!(p1, xi, maximum(yi)*Z, zi, color=:gray, markeralpha=markeralpha, markerstrokewidth=0)
		# actual points
		scatter3d!(p1, xi, yi, zi, m=:bluesreds_r, zcolor=zi, markerstrokewidth=0)
	end

	function createPlots(MODE, ii, scatcam; title="", ypl=nothing)
		xyzi, Nhead, ylabel = prepareDataNonlinBenefit(fname, MODE, ii)
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
				scatterAndShadow(Tractual, koratio, FLspec, scatcam, ylabel),
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
		# FIXME: low/high I (1 vs. 3 in the last element really seem to not make a big difference)
		# Mode 1: al;Qdt;phi;wingdens
		createPlots(1, ([2,3,4,5], [1,2], 2, [1,2,3]), (30,60); title="120deg")..., 
		createPlots(1, ([2,3,4,5], [1,2], 1, [1,2,3]), (30,60); title="90deg")...,
		# plot(createPlots(1, [4, 1, 2, 2])..., createPlots(1, [4, 1, 2, 3])..., title="Wing density (low, high)"),
	]
	# createPlots(1, [3, 1, 2, 3])
	# createPlots(2, [2,3,1,2])
end

"Normalize by lowering T1 corresponding to the maximum output disp"
function nlNormalizeByOutput(m, opt, param, σa0, T1scale; τ2ratio=2.0)
	if !isnothing(T1scale)
		param[2] *= T1scale
	else
		# Automatically try to match output. From transmission:
		# τfun = σa -> τ1*σa + τ2/3*σa^3
		# = τ1*σa + τ1*τ2ratio/3*σa^3 = τ1*σa*(1 + τ2ratio/3*σa^2)
		# So scale τ1 by 1/(1 + τ2ratio/3*σa^2)
		param[2] *= 1/(1 + τ2ratio/3*σa0^2)
	end
	param[5] = τ2ratio*param[2]
	return param
end

"""Generate plot like in [Jafferis (2016)] Fig. 4"""
function openLoopPlot(m, opt, param0, Vamps; save=false, NLT1scales=nothing, rightplot=false)
	function getResp(f, uamp, nlt, σa0=nothing, T1scale=nothing)
		param = copy(param0)
		if nlt
			param = nlNormalizeByOutput(m, opt, param, σa0, T1scale)
		end
		ts = createInitialTraj(m, opt, 0, f, [1e3, 1e2], param, 0; uampl=uamp, trajstats=true, thcoeff=0.1)
		# println("act disp=",ts[end])
		return ts
	end
	fs = 0.03:0.005:0.25
	mN_PER_V = 75/160

	p1 = plot(ylabel=rightplot ? "" : "Norm. stroke ampl [deg/V]", ylims=(0.3,0.8), legend=rightplot ? false : :topleft, title=rightplot ? "Low inertia" : "High inertia")
	p2 = plot(xlabel="Freq [kHz]", ylabel="Hinge ampl [deg]", legend=false, ylims=(0,100))
	p3 = plot(ylabel=rightplot ? "" : "Norm. act. disp [um/V]", legend=false, ylim=(1.0,2.4), xlabel="Freq [kHz]")
	# if rightplot
	# 	yaxis!(p1, false)
	# 	yaxis!(p3, false)
	# end
	function plotForTrans(nlt; T1scales=nothing)
		nltstr = nlt ? "N" : "L"
		actdisps = Dict{Float64, Float64}()
		
		for ii=1:length(Vamps)
			Vamp = Vamps[ii]
			T1scale = isnothing(T1scales) ? nothing : T1scales[ii]

			print("Openloop @ ", Vamp, "V ", nltstr)
			uamp = Vamp*mN_PER_V
			σa0 = nothing
			# println("HI", σa0)
			amps = hcat(getResp.(fs, uamp, nlt, σa0, T1scale)...)
			
			actdisps[Vamp] = maximum(amps[3,:])
			println(", act disp=", round(actdisps[Vamp], digits=3), "mm")
			amps[1:2,:] *= 180/pi # to degrees
			amps[1,:] /= (Vamp) # normalize
			amps[3,:] *= (1000/Vamp)
			amps[2,:] /= 2.0 # hinge ampl one direction
			# println(amps)
			plot!(p1, fs, amps[1,:], lw=2, label=string(nltstr, Vamp,"V"), ls=nlt ? :solid : :dash)
			plot!(p2, fs, amps[2,:], lw=2, label=string(nltstr, Vamp,"V"), ls=nlt ? :solid : :dash)
			plot!(p3, fs, amps[3,:], lw=2, label=string(nltstr, Vamp,"V"), ls=nlt ? :solid : :dash)
		end
		return actdisps
	end

	plotForTrans(false)
	plotForTrans(true; T1scales=NLT1scales)

	println("dens=", param0[3]/param0[6], ", koratio=", m.kbo[1]/(m.kbo[1] + m.ka/param0[2]^2))
	# plot(p1, #= p2, layout=(2,1),  =# size=(400, 300), dpi=200)
	return plot(p1, p3, layout=(2,1))
end

function openLoopPlotFinal(m, opt, param0)
	# mw=0.55 -> [0.975, 0.925], mw=0.7 -> [0.96, 0.9]
	param = copy(param0)
	param[3] = 0.55 #mw
	pl1 = openLoopPlot(m, opt, param, range(140, 180,length=2); NLT1scales=[0.975, 0.925])
	param[3] = 0.7 #mw
	pl2 = openLoopPlot(m, opt, param, range(140, 180,length=2); NLT1scales=[0.96, 0.9], rightplot=true)
	plot(pl1, pl2, size=(500, 500), dpi=150)
	gui()
	error("Open loop plot")
end
