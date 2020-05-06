using LinearAlgebra, DifferentialEquations, Plots, OSQP, SparseArrays, ForwardDiff

# TYPES ----------------------

abstract type ControlAffine end
struct CAVertical <: ControlAffine end
struct CAPlanar <: ControlAffine end

# ----------------------------

"As in the vertical case, v = Φ^(1/2)"
function aeroWrench(f, Φ2)
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

	Fz = kt * f^2 * Φ2 * cos(Ψ)*sin(Ψ)
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

# OSQP basic --------------

function qpSetupDense(n, m)
	model = OSQP.Model()
	# QP solution with only u
	# with a desired dydes, want 
	P = sparse(ones(n,n))
	q = zeros(n)
	l = fill(-Inf, m)
	u = fill(Inf, m)
	A = sparse(ones(m,n))
	OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, eps_rel=1e-2, eps_abs=1e-2, verbose=false)
	return model
end

# ----
function runSim(ca::ControlAffine, y0, tend, controller; simdt=0.1, udt=1)
	# Functions to evaluate the vector field
	lastUUpdate = -Inf
	Uts = []
	Us = []

	function vf(y, p, t)
		fy, gy = nonLinearDynamics(ca, y)

		if t > lastUUpdate + udt
			append!(Uts, t)
			push!(Us, controller(ca, t, udt, y, fy, gy))
			lastUUpdate = t
		end

		return fy + gy * Us[end]
	end
	# OL traj1
	teval = collect(0:simdt:tend) # [ms]
	prob = ODEProblem(vf, y0, (teval[1], teval[end]))
	sol = solve(prob, saveat=teval)
	yk = hcat(sol.u...) # ny,N
	uk = hcat(Us...)
	return sol.t, yk, Uts, uk
end

# test vertical 1D model ------------------

function ddq(m::CAVertical, y)
	z, dz, v = y
	mb = 100 #[mg]
	g = 9.81e-3 #[mN/mg]
	kaero = 5.0 # very approx: mN for v = 1
	return 1/mb * (kaero*v) - g
end

"ddz = v; dv = k(vdes - v)"
function nonLinearDynamics(m::CAVertical, y)
	k = 1.0 # first order vdot
	z, dz, v = y
	# control-affine cts
	dydt0 = [dz; ddq(m, y); -k * v]
	dydt1 = [0; 0; k] # * vdes (input)
	return dydt0, dydt1
end

cav = CAVertical()
y0 = zeros(3)
model = qpSetupDense(1, 1)

function cavController(ca, t, dt, y, fy, gy)
	zdotdes = 1 # zdotdes
	# discretized model ZOH. dydt = f(y) + g(y)v. y2 = y1 + dydt*dt = (y1 + dt * fy) + dt * dy * v
	fd = (y + dt * fy)
	gd = dt * gy
	# project by A1 into the only state needed by the objective
	A1 = ForwardDiff.gradient(yy -> ddq(ca, yy), y)
	a1 = ddq(ca, y) - dot(A1, y)
	ft = y[2] + a1 * dt + dot(A1, fd)
	gt = dot(A1, gd) # these are now just scalars
	P = [gt^2]
	q = [gt * (ft - zdotdes)]
	l = [0.0]
	u = [1.0]
	# update OSQP
	OSQP.update!(model, Px=P[:], q=q, l=l, u=u)
	# solve
	res = OSQP.solve!(model)
	return res.x[1] # since it is a scalar
end
tt, yy, tu, uu = runSim(cav, y0, 200, cavController; udt=2)

p2 = plot(tt, yy[2,:], lw=2, xlabel="t", label="dz", ylabel="dz", legend=false)
p3 = plot(tt, yy[3,:], lw=2, xlabel="t", label="v", ylabel="v")
plot!(p3, tu, uu[1,:], lw=2, xlabel="t", label="vdes")
plot(p2, p3)
gui()

## --- Planar model -----------

function ddq(m::CAPlanar, y)
	# TODO: freq in input
	fL = 0.15
	fR = 0.15
	# as in the vertical example, 
	vL = y[7]
	vR = y[8]

	return controlAffinePlanarDynamics(y[1:3], y[4:6], [fL, vL], [fR, vR])
end

"As in the vertical model, use v or Φ2"
function nonLinearDynamics(m::CAPlanar, y)
	# unpack
	vL = y[7]
	vR = y[8]
	# control-related
	kv = 1
	fns = deg2rad(0.5) # TODO: function of freq

	# Similar to vertical
	fy = [y[4:6];  ddq(m, y);  -kv * vL; -kv * vR]
	gy = [zeros(6, 2); diagm(0 => kv * [fns, fns])]

	return fy, gy
end

cap = CAPlanar()
y0 = zeros(8)#vcat(zeros(6),π/2,π/2)
ny = length(y0)
fy, gy = nonLinearDynamics(cap, y0)
# print(fieldnames(typeof(model.workspace)))

# model = qpSetupDense(2, 2)
# function capController(ca, t, dt, y, fy, gy)
# 	dqbdes = [0.0,1.0,0.0]
# 	# discretized model ZOH. dydt = f(y) + g(y)v. y2 = y1 + dydt*dt = (y1 + dt * fy) + dt * dy * v
# 	fd = (y + dt * fy)
# 	gd = dt * gy
# 	# project by A1 into the only states needed by the objective: dqb2 = dqb1 + dt * 
# 	A1 = [0, 1, dt]
# 	ft = dot(A1, fd)
# 	gt = dot(A1, gd) # these are now just scalars
# 	P = [gt^2]
# 	q = [gt' * (ft - dqbdes)]
# 	l = [0.0,0.0]
# 	u = [200.0,200.0]
# 	# update OSQP
# 	OSQP.update!(model, Px=P[:], q=q, l=l, u=u)
# 	# solve
# 	res = OSQP.solve!(model)
# 	return res.x[1] # since it is a scalar
# end
function capController(ca, t, dt, y, fy, gy)
	# vdes = [0,0.1,0]
	# return qpSolve(ca, model, fy, gy, vdes)
	return [150., 140.]
end
tt, yy, tu, uu = runSim(cap, y0, 100, capController; udt=2)
# vf(y0, [], 0)

# Plot
p1 = plot(yy[1,:], yy[2,:], lw=2, xlabel="y", ylabel="z", legend=false)
p2 = plot(tt, yy[4,:], lw=2, xlabel="t", label="dy")
plot!(p2, tt, yy[5,:], lw=2, xlabel="t", label="dz")
p3 = plot(tt, yy[3,:], lw=2, xlabel="t", label="phi")
# Plot inputs
pu1 = plot(tu, uu[1,:], lw=2, xlabel="t", label="VL")
plot!(pu1, tu, uu[2,:], lw=2, xlabel="t", label="VR")

plot(p1, p2, p3, pu1)
gui()
