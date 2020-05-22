using LinearAlgebra, DifferentialEquations, Plots, OSQP, SparseArrays, ForwardDiff

# TYPES ----------------------

abstract type ControlAffine end
struct CAVertical <: ControlAffine end
struct CARollPlane <: ControlAffine end

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

"N * y returns next pos, current vel"
function nextPos(nq, dt; inclVel=true)
	Ii = diagm(0 => ones(nq))
	return inclVel ? [Ii  dt*Ii; zeros(nq,nq)  Ii] : [Ii  dt*Ii]
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

function reactiveQPAffine(ca, y, dt)
	# current state
	fT0, gT0, fA0, gA0 = nonLinearDynamicsTAD(ca, y, dt)
	# Need to linearize at the next state for the low-res next dynamics
	nT = length(fT0)
	yA0 = y[nT+1:end]
	yT1 = fT0 + gT0 * yA0
	fT1, gT1 = nonLinearDynamicsTAD(ca, [yT1;yA0], dt)[1:2] # Only need fT,gT so does not matter what yA1 is
	# Return both pos, vel components
	ft = fT1 + gT1 * fA0
	gt = gT1 * gA0
	return ft, gt, yT1
end

## test vertical 1D model ------------------

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
	gT = reshape([0; a1], 2, 1)
	# Lower rows (anchor)
	fA = [-k * v]
	gA = reshape([k], 1, 1)
	return fT, gT, fA, gA
end

cav = CAVertical()
y0 = zeros(3)
model = qpSetupDense(1, 1)

function cavController(ca, t, dt, y)
	# current state
	# z0, dz0, yA0 = y
	ft, gt, yT1 = reactiveQPAffine(ca, y, dt)
	velControl = false
	if velControl
		ydesproj = [1.] # zdotdes
		# project to vel components
		ft = ft[2:2]
		gt = gt[2:2,:]
		wy = [1.]
	else
		# position control
		N = nextPos(1, dt)
		ft = N * ft
		gt = N * gt
		ydesproj = [10., 0.]
		wy = [1.,30.]
	end
	P = gt' * Diagonal(wy) * gt
	q = gt' * Diagonal(wy) * (ft - ydesproj)
	l = [0.0]
	u = [1.0]
	# update OSQP
	OSQP.update!(model, Px=P[:], q=q, l=l, u=u)
	# solve
	res = OSQP.solve!(model)
	
	return res.x # since it is a scalar
end
tt, yy, tu, uu = runSim(cav, y0, 200, cavController; udt=2)

p1 = plot(tt, yy[1,:], lw=2, xlabel="t", label="z", ylabel="z", legend=false)
p2 = plot(tt, yy[2,:], lw=2, xlabel="t", label="dz", ylabel="dz", legend=false)
p3 = plot(tt, yy[3,:], lw=2, xlabel="t", label="v", ylabel="v")
plot!(p3, tu, uu[1,:], lw=2, xlabel="t", label="vdes")
plot(p1, p2, p3)
gui()

## --- Planar model -----------

"When multiplied by Φ^2, gives the aero wrench"
function aeroWrenchAffine(f)
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

	Fz = kt * f^2 * cos(Ψ)*sin(Ψ)
	return [Fz, ycp*Fz]
end

"Returns a0, a1, s.t. ddq = a0 + a1 * u. Here u = Φ^2"
function caddq(m::CARollPlane, y)
	# TODO: freq in input
	fL = 0.15
	fR = 0.15

	# dynamics
	mb = 100 #[mg]
	ib = 3333 #[mg-mm^2]
	g = 9.81e-3 #[mN/mg]

	# dynamics stuff
	Mb = Diagonal([mb, mb, ib])
	h = [0; mb*g; 0]
	rot(x) = [cos(x) -sin(x); sin(x) cos(x)]
	# now should be multiplied by u = [ΦL^2; ΦR^2]
	totalWrenchAffine = hcat(aeroWrenchAffine(fL), Diagonal([1,-1]) * aeroWrenchAffine(fR)) # 2x2 matrix

	a0 = -Mb \ h # 3x1
	# println("HI ", y[3])
	a1 = Mb \ ([rot(y[3])[:,2] zeros(2,1); 0 1] * totalWrenchAffine) # this is now 3x2, should be multiplied by u = [ΦL^2; ΦR^2]
	# display(totalWrenchAffine)
	# display([rot(y[3])[:,1] zeros(2,1); 0 1])
	# display(a1)

	return a0, a1
end

"ddz = v; dv = k(vdes - v)"
function nonLinearDynamicsTA(m::CARollPlane, y)
	# unpack
	qb = y[1:3]
	dqb = y[4:6]
	yA = y[7:8]
	# control-related
	k = 1.0 # first order vdot
	fns = 1.#deg2rad(0.5)

	# control-affine cts
	a0, a1 = caddq(m, y)
	fT = [dqb; a0]
	gT = [zeros(3,2); a1]
	# Lower rows (anchor)
	fA = -k * yA
	gA = Diagonal(k * [fns, fns])
	return fT, gT, fA, gA
end

cap = CARollPlane()
y0 = zeros(8)#vcat(zeros(6),π/2,π/2)
ny = length(y0)
fy, gy = nonLinearDynamics(cap, y0)
# print(fieldnames(typeof(model.workspace)))

# REACTIVE CONTROLLER ---
model = qpSetupDense(2, 2)
function capController(ca, t, dt, y)
	velControl = false
	ft, gt, yT1 = reactiveQPAffine(ca, y, dt)
	Ax = Float64[1,0,0,1]

	if velControl
		ft = ft[4:6]
		gt = gt[4:6,:]
		# return [1.2,1]#deg2rad(0.5)*[150., 140.]
		wy = [0.,1.,1.]
		
		# # sine traj
		# y1des = sin(0.1*t)
		# dqbdes = [0.0,0.1,0.1*(y[1] - y1des)-1.0*y[3]]
		
		# circle traj
		y1des = 5*sin(0.01*t)
		ydesproj = [0.0,
			0.01*(5*(1+sin(0.01*t)) - y[2]),
			0.1*(y[1] - y1des)-1.0*y[3]]
	else
		# position control
		N = nextPos(3, dt)
		ft = N * ft
		gt = N * gt
		wy = [1.,1.,1.,50.,50.,20.]
		ydesproj = vcat([1.,10.,0.], zeros(3))
	end
	
	P = gt' * Diagonal(wy) * gt
	q = gt' * Diagonal(wy) * (ft - ydesproj)
	l = [0.0,0.0]
	u = [200.0,200.0]
	# update OSQP
	OSQP.update!(model, Px=P[triu!(trues(size(P)),0)], q=q, l=l, u=u, Ax=Ax)
	# solve
	res = OSQP.solve!(model)
	return res.x
end
# # MPC N=1 --------------------------
# # Increased size: now x=[uA0, uA1]
# nx = 4
# model = qpSetupDense(nx, nx)
# function capController(ca, t, dt, y)
# 	# return [1.2,1]#deg2rad(0.5)*[150., 140.]
# 	wy = [0.,1.,1.]
# 	Amat = Diagonal(ones(nx))
# 	Ax = Array(Amat[:])
	
# 	dqbdes = [0,0.1,0]#-1.0*y[3]]
# 	ft, gt, yT1 = reactiveQPAffine(ca, y, dt)
# 	# Need W2, i.e. first a predicted yT2
# 	yT2 = yT1 # TODO: other options
# 	W2 = caddq(ca, yT2)[2]
# 	# Also need affine form fA(y1)
# 	yA1 = y[7:8] # FIXME:
# 	fA1, gA1 = nonLinearDynamicsTAD(ca, [yT1; yA1], dt)[3:4]

# 	# Construct the QP
# 	A = ft + W2 * fA1
# 	B = hcat(gt, W2 * gA1)
# 	P = B' * Diagonal(wy) * B
# 	q = B' * Diagonal(wy) * (A - dqbdes)
# 	l = zeros(nx)
# 	u = 200.0 * ones(nx)
# 	# update OSQP
# 	OSQP.update!(model, Px=P[triu!(trues(size(P)),0)], q=q, l=l, u=u, Ax=Ax)
# 	# solve
# 	res = OSQP.solve!(model)
# 	return res.x[1:2]
# end

tt, yy, tu, uu = runSim(cap, y0, 500, capController; udt=2)
# vf(y0, [], 0)

# Plot
# p1 = plot(yy[1,:], yy[2,:], lw=2, xlabel="y", ylabel="z", legend=false)
p1 = plot(tt, yy[2:3,:]', lw=2, xlabel="t", ylabel="z", legend=false)
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
# savefig("test.png")
