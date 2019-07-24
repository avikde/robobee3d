
using ForwardDiff
include("Wing2DOF.jl")

println("hi")
y0 = [0.1,0.1,1,0]
u0 = [0.]
params0 = [0.05,1.0]
paeroFun(q::Vector) = Wing2DOF.aero([q;[0,0]], u0, params0)[1]
println("paero ", paeroFun(y0[1:2]))

# Try to use ForwardDiff
JaeroFun = y -> ForwardDiff.jacobian(paeroFun, y)
# println(JaeroFun(y0[1:2]))
println(Wing2DOF.aero(y0, u0, params0)[2])


println(Wing2DOF.dydt(y0, u0, params0))
