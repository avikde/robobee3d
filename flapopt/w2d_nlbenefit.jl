
NLBENEFIT_FNAME = "nonlin.zip"
function nonlinBenefit(ret, Tratios, minals)
	function maxu(τ21ratiolim, minal)
		rr = opt1(m, ret["traj"], ret["param"], 1, minal, τ21ratiolim)
		return [rr["u∞"]; rr["al"]; rr["FD∞"]; rr["param"]]
	end

	res = [[T;al;maxu(T,al)] for T in Tratios, al in minals]
	matwrite(NLBENEFIT_FNAME, Dict("res" => res); compress=true)
	return res
end

function plotNonlinBenefit(fname, ypl; s=100)
	results = matread(fname)["res"]
	
	# Row 1 has Tratio=0 (res[1,:])
	for r=size(results,1):-1:1
		for c=1:size(results,2)
			resT0 = results[1,c]
			# normalize
			results[r,c][3] /= resT0[3]
			results[r,c][5] /= resT0[5]
		end
	end
	xyzi = zeros(length(results[1,1]),0)
	for res in results
		xyzi = hcat(xyzi, res)
	end
	# lift to mg
	params = xyzi[6:end,:]

	xpl = [0,3]
	X = range(xpl[1], xpl[2], length=50)
	Y = range(ypl[1], ypl[2], length=50)

	function contourFromUnstructured(xi, yi, zi; title="")
		# Spline from unstructured data https://github.com/kbarbary/Dierckx.jl
		# println("Total points = ", length(xi))
		spl = Spline2D(xi, yi, zi; s=s)
		ff(x,y) = spl(x,y)
		return contour(X, Y, ff, 
			titlefontsize=10, grid=false, lw=2, c=:bluesreds, 
			xlabel="T ratio", ylabel="FL [mg]", title=title,
			xlims=xpl, ylims=ypl)
	end
    
	return [
		# scatter(xyzi[1,:], xyzi[4,:]),
		# scatter3d(xyzi[1,:], xyzi[4,:], xyzi[3,:]),
		contourFromUnstructured(xyzi[1,:], xyzi[4,:], clamp.(xyzi[3,:], 0.0, 1.0); title="Nonlinear transmission benefit [ ]"),# opt errors cause > 1 which doesn't make sense
		contourFromUnstructured(xyzi[1,:], xyzi[4,:], params[2,:]; title="T1 [rad/mm]"),
		contourFromUnstructured(xyzi[1,:], xyzi[4,:], params[6,:]; title="Aw [mm^2]"),
		contourFromUnstructured(xyzi[1,:], xyzi[4,:], 1000.0 ./(N*params[7,:]); title="Freq [Hz]")
	]
end
