
using DifferentialEquations
using Plots; gr()

import controlutils
cu = controlutils

mutable struct MassSpringDamperModel <: controlutils.Model
    # Actuator params
    ma::Float64 # [mg]; effective
    ka::Float64 # [mN/mm]
    mo::Float64 # [mg]
    umax::Float64
end

function cu.dims(m::MassSpringDamperModel)::Tuple{Int, Int}
    return 2, 1
end

function cu.pdims(m::MassSpringDamperModel)::Int
    return 1
end

function cu.plimits(m::MassSpringDamperModel)
    return [0.1, 1.0]
end

function cu.limits(m::MassSpringDamperModel)::Tuple{Vector, Vector, Vector, Vector}
    # This is based on observing the OL trajectory. See note on units above.
    umax = [m.umax] # [mN]
    umin = -umax
    xmax = [10,1000] # [mm, mm/ms]
    xmin = -xmax
    return umin, umax, xmin, xmax
end

function cu.limitsTimestep(m::MassSpringDamperModel)::Tuple{Float64, Float64}
	return 0.01, 0.07
end

"Continuous dynamics second order model"
function cu.dydt(m::MassSpringDamperModel, y::AbstractArray, u::AbstractArray, param::Vector)::AbstractArray
    # unpack
    τ1, ko, bo, τ2 = param
    ya, T, τfun, τifun = cu.transmission(m, y, param; o2a=true)
    σo, dσo = y
    ddy = 1.0/(m.mo + m.ma/T^2) * (-(ko*σo + m.ka/T*τifun(σo)) - (bo)*dσo + u[1]/T)
    # return ddq
    return [y[2], ddy]
end

# "Objective to minimize"
# function cu.robj(m::MassSpringDamperModel, opt::cu.OptOptions, traj::AbstractArray, params::AbstractArray)::AbstractArray
#     ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    
# 	yk = k -> @view traj[liy[:,k]]
# 	uk = k -> @view traj[liu[:,k]]
#     kvec = collect(1:N+1)
#     rtrajErr = first.(yk.(kvec)) - σdes.(N, kvec)
#     # only traj cost
#     return rtrajErr
#     # # also add input cost?
#     # ru = 0.1 * uk.(collect(1:N))
#     # return [rtrajErr;ru]
# end

# -------------------------------------------------------------------------------


function createOLTraj(m::MassSpringDamperModel, opt::cu.OptOptions, traj::AbstractArray, params::AbstractArray)::AbstractArray
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    
    traj1 = copy(traj)
    yk = k -> @view traj1[liy[:,k]]
	uk = k -> @view traj1[liu[:,k]]
    for k=1:N
        traj1[liy[:,k+1]] = cu.ddynamics(m, yk(k), [0], params, δt)
        traj1[liu[:,k]] = [0]
    end
    return traj1
end

function refTraj(m::MassSpringDamperModel, freq)
    σmax = cu.limits(m)[end][1]
    # had trouble with
    post = t -> σmax * cos(freq * 2 * π * t)
    velt = t -> -σmax * (freq * 2 * π) * sin(freq * 2 * π * t)
    acct = t -> -σmax * (freq * 2 * π)^2 * cos(freq * 2 * π * t)
    return post, velt, acct
end

"""
freq [kHz]; posGains [mN/mm, mN/(mm-ms)]; [mm, 1]
Example: trajt, traj0 = Wing2DOF.createInitialTraj(0.15, [1e3, 1e2], params0)
"""
function createInitialTraj(m::MassSpringDamperModel, opt::cu.OptOptions, N::Int, freq::Real)
    # Create an output traj (came from a template or something)
    post, velt = refTraj(m, freq)[1:2]
    trajt = 0:opt.fixedδt:(N)*opt.fixedδt
    # dirtran form {x1,..,x(N+1),u1,...,u(N),δt}
    dirtrany = reshape(hcat(post.(trajt), velt.(trajt))', :)
    traj = [dirtrany; zeros(N)]
    return trajt, traj, trajt
end

function plotTrajs(m::MassSpringDamperModel, opt::cu.OptOptions, t, params, trajs; ulim=nothing, fdes=0.1)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, trajs[1])
    Ny = (N+1)*ny
    # If plot is given a matrix each column becomes a different line
    σt = plot(t, hcat([traj[@view liy[1,:]] for traj in trajs]...), linewidth=2, ylabel="stroke [mm]", title="δt=$(round(δt; sigdigits=4))ms")
    dσt = plot(t, hcat([traj[@view liy[2,:]] for traj in trajs]...), linewidth=2, ylabel="stroke vel [mm/ms]", legend=false)
    
    ut = plot(t, hcat([[traj[@view liu[1,:]];NaN] for traj in trajs]...), linewidth=2, legend=false, ylabel="stroke force [mN]")
    if !isnothing(ulim)
        ylims!(ut, (-ulim, ulim))
    end
    
    # actuator displacement
    function actdisp(traj, param)
        τ1, ko, bo, τ2 = param
        y2, T, τfun, τifun = cu.transmission(m, traj, param; o2a=true)
        return τifun.(traj[@view liy[1,:]])
    end
    σat = plot(t, hcat([actdisp(trajs[i], params[i]) for i=1:length(trajs)]...), linewidth=2, ylabel="act stroke [mm]")

    # also plot the des pos and vel to make sure the initial traj is "OK"
    post, velt = refTraj(m, fdes)
    plot!(σt, t, post.(t), color=:black, linestyle=:dash)
    plot!(dσt, t, velt.(t), color=:black, linestyle=:dash)
    
    # Combine the subplots
	return (σt, dσt, σat, ut)
end

# ---------------------- Param opt -----------------------------

function cu.paramLumped(m::MassSpringDamperModel, param::AbstractArray)
    τ1, ko, bo, τ2 = param
    return [1, ko, bo], [τ1, τ2]
end

function cu.transmission(m::MassSpringDamperModel, y::AbstractArray, _param::Vector; o2a=false)
    τ1, ko, bo, τ2 = _param
    τfun = σa -> τ1*σa + τ2/3*σa^3
    # Series from Mathematica
    τifun = σo -> σo/τ1 - τ2*σo^3/(3*τ1^4) # + O[σo^4]
    if !o2a
        Dτfun = σa -> τ1 + τ2*σa^2
        T = Dτfun(y[1]) # "gear ratio"
        y2 = [τfun(y[1]), T*y[2]] # [mm, mm/ms]
    else
        σo = y[1]
        T = τ1 + σo^2*τ2/τ1^2
        y2 = [τifun(σo), y[2]/T] # [mm, mm/ms]
    end
    return y2, T, τfun, τifun
end


function cu.paramAffine(m::MassSpringDamperModel, opt::cu.OptOptions, traj::AbstractArray, param::AbstractArray, POPTS::cu.ParamOptOpts, scaleTraj=1.0; debugComponents=false)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)

    yo = k -> @view traj[liy[:,k]] # traj is in output coords already
    uk = k -> @view traj[liu[:,k]]

    # Param stuff
    τ1, ko, bo, τ2 = param

    # THESE FUNCTIONS USE OUTPUT COORDS -------------
    function HMqTo(ypos, yvel)
        σo, σ̇odum = ypos * scaleTraj
        σodum, σ̇o = yvel * scaleTraj
        return [σ̇o*m.mo   0   0   0   0]
    end
    function HMqTa(ypos, yvel)
        σo, σ̇odum = ypos * scaleTraj
        σodum, σ̇o = yvel * scaleTraj
        return [0   0   0   σ̇o*m.ma   σ̇o*m.ma*(-σo^2)]
    end
    HMqT(ypos, yvel) = HMqTo(ypos, yvel) + HMqTa(ypos, yvel)
    Hio = (y, ynext) -> HMqTo(y, ynext) - HMqTo(y, y)
    Hia = (y, ynext) -> HMqTa(y, ynext) - HMqTa(y, y)

    function Hstiffo(y)
        σo, σ̇o = y * scaleTraj
        return δt*[0   σo   0   0   0]
    end
    function Hstiffa(y)
        σo, σ̇o = y * scaleTraj
        return δt*[0   0   0   m.ka*σo   m.ka*(-σo^3/3)]
    end
    function Hdamp(y)
        σo, σ̇o = y * scaleTraj
        return δt*[0   0   σ̇o   0    0]
    end
    HCgJTδt(y) = Hstiffo(y) + Hstiffa(y) + Hdamp(y)

    if debugComponents
        function Hτ(Hh, y)
            σo = y[1]
            return hcat(Hh[:,1:end-2], Hh[:,1:end-2]*σo^2, Hh[:,end-1:end])
        end
        return yo, Hτ, Hio, Hia, Hstiffo, Hstiffa, Hdamp
    end
    # ----------------

    # 1st order integration
    Htil = (y, ynext) -> HMqT(y, ynext) - HMqT(y, y) + HCgJTδt(y)
    
    # Functions to output
    "Takes in a Δy in output coords"
    function Hk(k, Δyk, Δykp1)
        Hh = Htil(yo(k) + Δyk, yo(k+1) + Δykp1)
        # With new nonlinear transmission need to break apart H
        σo = (yo(k) + Δyk)[1]
        return hcat(Hh[:,1:end-2], Hh[:,1:end-2]*σo^2, Hh[:,end-1:end])
    end
    
    # For a traj, H(yk, ykp1, Fk) * pb = B uk for each k
    B = reshape([1.0], 1, 1)

    return Hk, yo, uk, B, N
end

