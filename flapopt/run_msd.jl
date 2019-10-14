
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
using SparseArrays, OSQP, LinearAlgebra # temp
import controlutils
cu = controlutils
includet("MassSpringDamper.jl")

# create an instance
m = MassSpringDamperModel(6, # ma
	150, # ka
	500, # mo
	100) # umax
ny, nu = cu.dims(m)
opt = cu.OptOptions(false, 0.2, 1, :none, 1e-8, false)
N = 34
param0 = [10.0, # T
100.0, # ko
10.0] # bo

fdes = 0.15
# FIXME: template
trajt, traj0, trajt = createInitialTraj(m, opt, N, fdes, [1e3, 1e2], param0)

# Make traj satisfy dyn constraint with these params?
traj0 = cu.fixTrajWithDynConst(m, opt, traj0, param0)


R_WTS = (zeros(2,2), 0, 1.0*I)#diagm(0=>[0.1,100]))
Tmin = 10.0 # FIXME: calculate
plimsL = [Tmin, 0.1, 0.1]
plimsU = [100.0, 100.0, 100.0]

param1, _, traj1, unactErr = cu.optAffine(m, opt, traj0, param0, 1, R_WTS, 0.1, plimsL, plimsU; Fext_pdep=true, test=true, testTrajReconstruction=false, print_level=1, max_iter=100)
# display(param1')


# # Optimize params directly for traj
# Hh = [Hdes(m, fdes, t) for t in trajt1]
# np = 3
# P1 = zeros(np,np)
# for Hk in Hh
# 	P1 .= P1 + Hk * Hk'
# end
# mo = OSQP.Model()
# OSQP.setup!(mo; P=sparse(P1), q=zeros(np), A=sparse(1:np, 1:np, ones(np)), l=[0; 0; m.mb], u=[Inf; m.bσ; m.mb])
# res = OSQP.solve!(mo)
# resf = sqrt(res.x[1]/m.mb)/(2*π)
# println("Res freq = ", resf)

# # Plot the quadratic
# size = 100
# x = range(0, stop=20, length=size)
# y = range(0, stop=10, length=size)
# ff(x,y) = [x;y;m.mb]' * P1 * [x;y;m.mb]
# pl1 = contour(x, y, ff)
# vline!(pl1, [res.x[1]])
