using Plots

# # Test feasibility
# function paramTest(p, paramConstraint)
# 	xtest = [p; zeros((N+1)*ny)]
# 	g = paramConstraint(xtest)
# 	gunact = g[1:(N+1)]
# 	grest = g[(N+1)+1:end]
# 	# have [gpolycon (2); gtransmission (1); ginfnorm (2*N)]
# 	println("Feas: should be nonpos: ", maximum(grest), "; unact: ", maximum(abs.(gunact)) ,"; transmission: ", g[3])
# 	unew = cu.getTrajU(m, opt, traj1, p, POPTS)
# 	println("Obj: ", paramObj(xtest))
# end

"See https://github.com/avikde/robobee3d/pull/136.
Also try to see how the traj was changed https://github.com/avikde/robobee3d/issues/142"
function debugΔYEffect(N, ny, rr, rr0=nothing)
	pt, Hk, B, Js, actVec = rr["eval_f"](rr["x"]; debug=true)
	println("Js ", Js)
	Δy = rr["x"][length(rr["param"])+1:end]
	Δyk(k) = Δy[(k-1)*ny+1:(k)*ny]
	Δy0 = zeros(ny)
	unew = vcat([B' * Hk(k,Δyk(max(k-1,1)),Δyk(k),Δyk(k+1))[1] * pt for k=1:N]...)
	unew0 = vcat([B' * Hk(k,Δy0,Δy0,Δy0)[1] * pt for k=1:N]...)
	p1 = plot([rr["traj"][(N+1)*ny+1:end]  actVec[:,1]  unew], lw=2, ls=[:solid :solid :dash])
	plot!(p1, unew0, lw=2, ls=:dash)

	#Infeasibility
	infeas = zeros(rr["nc"])
	rr["eval_g"](rr["x"], infeas)
	np = length(rr["param"])

	#Hinge deltay
	function ΔyOutputChange(i)
		pp = plot(rr["traj"][i:ny:(N+1)*ny], lw=2, ylabel=(i==1 ? "Stroke [rad]" : "Pitch [rad]"), label="new", legend=(i == 2 ? :bottomright : false))
		if !isnothing(rr0)
			plot!(pp, rr0["traj"][i:ny:(N+1)*ny] * (i==1 ? rr0["param"][end]/rr["param"][end] : 1), lw=2, label="orig")
		end
		return pp
	end

	D = cu.trajDiffMat(N, ny)

	return (p1, 
		plot([plot(Δy[i:ny:(N+1)*ny], legend=false, lw=2, ylabel=string(i)) for i=1:ny]...), 
		plot(D * rr["x"][np+1:end],lw=2,ylabel="diff(Δy)", legend=false), 
		plot(infeas[1:N],lw=2,ylabel="unact constraint", legend=false), 
		plot(infeas[N+1:N+np], ylims=(-0.1,0.1),lw=2,ylabel="polytope constraint", legend=false),
		plot(ΔyOutputChange(1), ΔyOutputChange(2)#= , layout=(2,1) =#)
		)
end

function debug4()
	function d1(minal, phi)
		ret2 = @time opt1(m, ret1["traj"], ret1["param"], 1, minal; Φ=phi)#, print_level=3)
		plot(debugComponentsPlot(m, opt, POPTS, ret2)..., debugDeltaYEffect(N, ny, ret2)..., size=(1600,800))
		savefig(string("debug4_", minal, "_", phi, ".png"))
	end
	d1(180,90)
	d1(180,120)
	d1(400,90)
	d1(400,120)
	error("debug4")
end

"https://github.com/avikde/robobee3d/pull/137"
function debugGradient(rr)
	x = copy(rr["x"])
	x1s = 5.0:40
	fx1(x1, a=true) = rr["eval_f"]([x1; x[2:end]]; auto=a)
	function ddx1(x1, a=true)
		grad_f = similar(x)
		rr["eval_grad_f"]([x1; x[2:end]], grad_f; auto=a)
		return grad_f[1]
	end
	plot(
		plot([fx1.(x1s) fx1.(x1s, false)],lw=2,ls=[:solid :dash]), 
		plot([ddx1.(x1s)  ddx1.(x1s, false)],lw=2,ls=[:solid :dash])
	)
	gui()
	error("h")
end

# function compareTrajToDAQ(m::Wing2DOFModel, opt::cu.OptOptions, t::Vector, param, traj, lift, drag)
# 	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
# 	Ny = (N+1)*ny
#     cbar2, τ1, mwing, wΨ, τ2, Aw, dt = param
#     # Plot of aero forces at each instant
#     function aeroPlotVec(_traj::Vector, _param, i)
#         Faerok = k -> w2daero(m, _traj[@view liy[:,k]], _param)[end]
#         Faeros = hcat([Faerok(k) for k=1:N]...)
#         return [Faeros[i,:]' NaN]'
#     end
#     liftp = plot(t, lift, marker=:auto, legend=false, label="daq", ylabel="lift [mN]")
#     plot!(liftp, t, aeroPlotVec(traj, param, 2), marker=:auto, legend=false, label="pred")
#     dragp = plot(t, drag, marker=:auto, label="daq", ylabel="drag [mN]")
#     plot!(dragp, t, aeroPlotVec(traj, param, 1), marker=:auto, label="pred")

#     # Get the basic kinematics
#     σt, Ψt, ut, _ = plotTrajs(m, opt, [param], [traj])

#     # Combine the subplots
# 	return (σt, Ψt, ut, liftp, dragp)
# end
