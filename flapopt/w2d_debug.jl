
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

"See https://github.com/avikde/robobee3d/pull/136"
function debugDeltaYEffect(rr)
	pt, Hk, B, Js, actVec = rr["eval_f"](rr["x"]; debug=true)
	println("Js ", Js)
	dy = rr["x"][length(param0)+1:end]
	dely(k) = dy[(k-1)*ny+1:(k)*ny]
	dely0 = zeros(ny)
	unew = vcat([B' * Hk(k,dely(k),dely(k+1))[1] * pt for k=1:N]...)
	unew0 = vcat([B' * Hk(k,dely0,dely0)[1] * pt for k=1:N]...)
	p1 = plot([rr["traj"][(N+1)*ny+1:end]  actVec[:,1]  unew], lw=2, ls=[:solid :solid :dash])
	plot!(p1, unew0, lw=2, ls=:dash)
	infeas = ret2["eval_g"](ret2["x"])
	return (p1, 
		plot(
			plot(dy[1:ny:(N+1)*ny]),
			plot(dy[2:ny:(N+1)*ny]),
			plot(dy[3:ny:(N+1)*ny]),
			plot(dy[4:ny:(N+1)*ny])
		), 
		plot(infeas[1:N],lw=2,ylabel="unact constraint"), 
		plot(infeas[N+1:end], ylims=(-0.1,0.1),lw=2,ylabel="polytope constraint")
		)
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
