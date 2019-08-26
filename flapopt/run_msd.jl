
# This is a "script" to test/run the functions from
# push!(LOAD_PATH, pwd()) # Only needs to be run once

using BenchmarkTools
using Revise # while developing
import controlutils
cu = controlutils
includet("MassSpringDamper.jl")

# create an instance
m = MassSpringDamperModel()
ny, nu = cu.dims(m)
opt = cu.OptOptions(false, 0.2, 1, :symmetric, 1e-8, false)
# opt = cu.OptOptions(false, 0.2, 1, cu.SYMMETRIC, 1e-8, false)
N = opt.boundaryConstraint == :symmetric ? 17 : 34
param0 = [1.0] # k

trajt, traj0 = createInitialTraj(m, opt, N, 0.15, [1e3, 1e2], param0)


# IPOPT
εs = [1, 0.005, 0.001] # IC, dyn, symm
prob = cu.ipoptsolve(m, opt, traj0, param0, εs, :traj; print_level=1, nlp_scaling_method="none")
traj1 = prob.x
# with my modification to Ha/Coros g-preferred param opt
δx = cu.paramδx(m, opt, traj1, param0, prob.mult_x_L, prob.mult_x_U)
param1 = cu.paramopt(nothing, m, opt, traj1, param0, δx, εs; step=1e2)
# param1 = cu.paramoptJ(m, opt, traj1, param0, εs; step=0.01)
prob = cu.ipoptsolve(m, opt, traj1, param1, εs, :traj; print_level=1, nlp_scaling_method="none")
traj2 = prob.x
# with my modification to Ha/Coros g-preferred param opt
δx = cu.paramδx(m, opt, traj2, param1, prob.mult_x_L, prob.mult_x_U)
param2 = cu.paramopt(nothing, m, opt, traj2, param1, δx, εs; step=1e2)
# param1 = cu.paramoptJ(m, opt, traj1, param0, εs; step=0.01)
prob = cu.ipoptsolve(m, opt, traj2, param2, εs, :traj; print_level=1, nlp_scaling_method="none")
traj3 = prob.x
# with my modification to Ha/Coros g-preferred param opt
δx = cu.paramδx(m, opt, traj3, param2, prob.mult_x_L, prob.mult_x_U)
param3 = cu.paramopt(nothing, m, opt, traj3, param2, δx, εs; step=1e2)
# param1 = cu.paramoptJ(m, opt, traj1, param0, εs; step=0.01)
prob = cu.ipoptsolve(m, opt, traj3, param3, εs, :traj; print_level=1, nlp_scaling_method="none")
traj4 = prob.x

# Results ---

trajs = [traj0, traj1, traj2, traj3, traj4]
params = [param0, param0, param1, param2, param3]

println("Objectives: ", [(_ro = cu.robj(m, opt, trajs[i], params[i]); _ro'*_ro) for i = 1:length(trajs)])
println("Params: ", params)

pl1 = plotTrajs(m, opt, trajt, params, trajs)
pl2 = cu.visualizeConstraintViolations(m, opt, params, trajs)

l = @layout [grid(2,1) a]
plot(pl1..., pl2, layout=l, size=(900,400))
