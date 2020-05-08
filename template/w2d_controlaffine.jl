using LinearAlgebra, DifferentialEquations, Plots, OSQP, SparseArrays, ForwardDiff

# TYPES ----------------------

abstract type ControlAffine end
struct CAVertical <: ControlAffine end
struct CAPlanar <: ControlAffine end

"Cascaded A->T system dynamics. dyT = fT + gT * yA; dyA = fA + gA * u. Returns fy, gy s.t. dy = fy + gy * u"
function nonLinearDynamics(m::ControlAffine, y)
	fT, gT, fA, gA = nonLinearDynamicsTA(m, y)
	nT = length(fT)
	nU = size(gA, 2)
	yA = y[nT+1:end]
	return [fT + gT * yA; fA], [zeros(nT, nU); gA]
end

"ZOH version"
function nonLinearDynamicsTAD(m::ControlAffine, y, dt)
	fT, gT, fA, gA = nonLinearDynamicsTA(m, y)
	# y1 = y0 + dt * dyT = y0 + dt * (fT + gT * yA) = (y0 + dt * fT) + (dt * gt) * yA
	nT = length(fT)
	yT0 = y[1:nT]
	yA0 = y[nT+1:end]
	# yA1 = yA0 + dt * dyA = (yA0 + dt * fA) + (dt * gA) * uA
	return yT0 + dt * fT, dt * gT, yA0 + dt * fA, dt * gA
end

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
			push!(Us, controller(ca, t, udt, y))
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

"Returns a0, a1, s.t. ddq = a0 + a1 * u"
function caddq(m::CAVertical, y)
	z, dz, v = y
	mb = 100 #[mg]
	g = 9.81e-3 #[mN/mg]
	kaero = 5.0 # very approx: mN for v = 1
	return -g, 1/mb * kaero
end

"ddz = v; dv = k(vdes - v)"
function nonLinearDynamicsTA(m::CAVertical, y)
	k = 1.0 # first order vdot
	z, dz, v = y
	# control-affine cts
	a0, a1 = caddq(m, y)
	fT = [dz; a0]
	gT = [0; a1]
	# Lower rows (anchor)
	fA = [-k * v]
	gA = [k]
	return fT, gT, fA, gA
end

cav = CAVertical()
y0 = zeros(3)
model = qpSetupDense(1, 1)

function cavController(ca, t, dt, y)
	zdotdes = [1.] # zdotdes
	# current state
	z0, dz0, yA0 = y
	fT0, gT0, fA0, gA0 = nonLinearDynamicsTAD(ca, y, dt)
	# Need to linearize at the next state for the low-res next dynamics
	yT1 = fT0 + gT0 * yA0
	fT1, gT1 = nonLinearDynamicsTAD(ca, [yT1;0], dt)[1:2] # Only need fT,gT so does not matter what yA1 is
	# yT2depu = gT1 * (fA0 + gA0 * u)
	# project by A1 into the only state needed by the objective
	ft = gT1[2] * fA0
	gt = gT1[2] * gA0
	P = [gt' * gt]
	q = [gt' * (ft - zdotdes)]
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
	fns = 1.0#deg2rad(0.5) # TODO: function of freq

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

model = qpSetupDense(2, 2)
function capController(ca, t, dt, y, fy, gy)
	dqbdes = [0.0,1.0,0.0]
	# discretized model ZOH. dydt = f(y) + g(y)v. y2 = y1 + dydt*dt = (y1 + dt * fy) + dt * dy * v
	fd = (y + dt * fy)
	gd = dt * gy
	# project by A1 into the only state needed by the objective
	A1 = ForwardDiff.jacobian(yy -> ddq(ca, yy), y)
	a1 = ddq(ca, y) - A1 * y
	ft = y[4:6] + (a1 + A1 * fd) * dt
	gt = dt * A1 * gd # these are now just scalars
	P = gt' * gt
	q = gt' * (ft - dqbdes)
	l = [0.0,0.0]
	u = [200.0,200.0]
	# update OSQP
	OSQP.update!(model, Px=P[:], q=q, l=l, u=u)
	# solve
	res = OSQP.solve!(model)
	return res.x # since it is a scalar
end
# function capController(ca, t, dt, y, fy, gy)
# 	# vdes = [0,0.1,0]
# 	# return qpSolve(ca, model, fy, gy, vdes)
# 	return [150., 140.]
# end
tt, yy, tu, uu = runSim(cap, y0, 20, capController; udt=2)
# vf(y0, [], 0)

# Plot
p1 = plot(yy[1,:], yy[2,:], lw=2, xlabel="y", ylabel="z", legend=false)
p2 = plot(tt, yy[4,:], lw=2, xlabel="t", label="dy")
plot!(p2, tt, yy[5,:], lw=2, xlabel="t", label="dz")
p3 = plot(tt, yy[3,:], lw=2, xlabel="t", label="phi")
# Plot inputs
pu1 = plot(tu, uu[1,:], lw=2, xlabel="t", ls=:dash, label="VLdes")
plot!(pu1, tt, yy[7,:], lw=2, xlabel="t", label="VL")
plot!(pu1, tu, uu[2,:], lw=2, xlabel="t", ls=:dash, label="VRdes")
plot!(pu1, tt, yy[8,:], lw=2, xlabel="t", label="VR")

plot(p1, p2, p3, pu1)
gui()
