using LinearAlgebra, DifferentialEquations, Plots, OSQP, SparseArrays

# control-affine
function aeroWrench(f, Φ)
	# params?
    CLmax = 1.8
    CDmax = 3.4
    CD0 = 0.4
	ρ = 1.225e-3 # [mg/(mm^3)]
	Aw = 54.4 #mm^2
	cbar = 3.2 #mm
	r2h = 0.551

	# calculated
	AR = Aw/cbar^2
	k = 1/2 * ρ * Aw^2 * AR * r2h^2
	kt = k * CLmax * π^2
	Lw = Aw/cbar
	ycp = 0.5 * Lw
	
	Ψ = 1.0 # TODO:

	Fz = kt * f^2 * Φ^2 * cos(Ψ)*sin(Ψ)
	return [Fz, ycp*Fz]
end

function controlAffinePlanarDynamics(qb, dqb, xwL, xwR)
	mb = 100 #[mg]
	ib = 3333 #[mg-mm^2]
	g = 9.81e-3 #[mN/mg]

	# dynamics stuff
	Mb = diagm(0 => [mb, mb, ib])
	h = [0; mb*g; 0]
	rot(x) = [cos(x) -sin(x); sin(x) cos(x)]
	totalWrench = aeroWrench(xwL...) + diagm(0=>[1,-1]) * aeroWrench(xwR...)
	Fz, Rx = totalWrench
	return Mb \ (-h + [rot(qb[3]) * [0;Fz]; Rx])
end

function controlAffinePlanar(y)
	# unpack
	qb = y[1:3]
	dqb = y[4:6]
	fL = y[7] # freq
	ΦL = y[8]
	fR = y[9] # freq
	ΦR = y[10]

	# control-related
	kf = 1
	kv = 1
	fns(f) = deg2rad(0.5) # TODO:

	ddq = controlAffinePlanarDynamics(qb, dqb, [fL, ΦL], [fR, ΦR])

	fy = [dqb; 
		ddq;
		-kf * fL;
		-kv * ΦL;
		-kf * fR;
		-kv * ΦR]
	
	gy = [zeros(6, 4);
		diagm(0 => [kf, kv * fns(fL), kf, kv * fns(fR)])]

	return fy, gy
end

# --- OSQP

model = OSQP.Model()
# QP solution with only u
n = 4
m = 4
P = sparse(ones(n,n))
q = zeros(n)
l = fill(-Inf, m)
u = fill(Inf, m)
A = sparse(ones(m,n))
OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, eps_rel=1e-2, eps_abs=1e-2)
res = OSQP.solve!(model)

# ----

# test
y0 = vcat(zeros(6),0.15,π/2,0.15,π/2)
fy, gy = controlAffinePlanar(y0)

function runSim(y0, tend; simdt=0.02)
	function vf(y, p, t)
		fy, gy = controlAffinePlanar(y)
		return fy + gy * [0.15, π/2, 0.15, π/2]
	end
	# OL traj1
	teval = collect(0:simdt:tend) # [ms]
	prob = ODEProblem(vf, y0, (teval[1], teval[end]))
	sol = solve(prob, saveat=teval)
	return sol.t, hcat(sol.u...)
end

tt, yy = runSim(y0, 1)
# vf(y0, [], 0)

# Plot
p1 = plot(yy[1,:], yy[2,:], lw=2, xlabel="y", ylabel="z", legend=false)
p2 = plot(tt, yy[1,:], lw=2, xlabel="t", label="y")
plot!(p2, tt, yy[2,:], lw=2, xlabel="t", label="z")
p3 = plot(tt, yy[3,:], lw=2, xlabel="t", label="phi")

plot(p1, p2, p3)
gui()
