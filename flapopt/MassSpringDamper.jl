
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
    τ1, ko, bo, τ2, dt = param
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

function plotTrajs(m::MassSpringDamperModel, opt::cu.OptOptions, t, params, trajs; ulim=nothing, fdes=0.1, refScale=1.0)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, trajs[1])
    Ny = (N+1)*ny
    Nt = length(trajs)
    
    # Empty versions of all the subplots
    σt = plot(ylabel="stroke [mm]")
    dσt = plot(ylabel="stroke vel [mm/ms]")
    σat = plot(ylabel="stroke act [mm]")
    ut = plot(ylabel="stroke force [mN]")
    if !isnothing(ulim)
        ylims!(ut, (-ulim, ulim))
    end
    
    # actuator displacement
    function actdisp(traj, param)
        τ1, ko, bo, τ2, dt = param
        y2, T, τfun, τifun = cu.transmission(m, traj, param; o2a=true)
        return τifun.(traj[@view liy[1,:]])
    end

    for i=1:Nt
        traj, param = trajs[i], params[i]
        dt = param[end]
        t = 0:dt:(N)*dt
        plot!(σt, t, traj[1:ny:(N+1)*ny], linewidth=2)
        plot!(dσt, t, traj[2:ny:(N+1)*ny], linewidth=2)
        plot!(σat, t, actdisp(traj, param), linewidth=2)
        plot!(ut, t, [traj[@view liu[1,:]];NaN], linewidth=2)
        # also plot the des pos and vel to make sure the initial traj is "OK"
        post, velt = refTraj(m, fdes)
        plot!(σt, t, refScale*post.(t), color=:blue, linestyle=:dash, lw=2)
        # FIXME: this factor seems like it is related to the freq. initial dt=0.1; if the dt is changed so that new dt=0.3; this needs to be newdt/olddt. Why??
        plot!(dσt, t, refScale*velt.(t), color=:blue, linestyle=:dash, lw=2)
    end

    # Combine the subplots
	return (σt, dσt, σat, ut)
end

# ---------------------- Param opt -----------------------------

function cu.paramLumped(m::MassSpringDamperModel, param::AbstractArray)
    τ1, ko, bo, τ2, dt = param
    return [1, ko, bo], τ1, τ2, dt
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

    yo = k -> traj[liy[:,k]]*scaleTraj # traj is in output coords already
    uk = k -> @view traj[liu[:,k]]

    # Param stuff
    τ1, ko, bo, τ2, dtold = param

    # THESE FUNCTIONS USE OUTPUT COORDS -------------
    function HMqTo(ypos, yvel)
        σo = ypos[1] * scaleTraj
        σ̇o = yvel[2] * scaleTraj * dtold
        return [σ̇o*m.mo   0   0   0   0]
    end
    function HMqTa(ypos, yvel)
        σo = ypos[1] * scaleTraj
        σ̇o = yvel[2] * scaleTraj * dtold
        return [0   0   0   σ̇o*m.ma   σ̇o*m.ma*(-σo^2)]
    end
    HMqT(ypos, yvel) = HMqTo(ypos, yvel) + HMqTa(ypos, yvel)
    Hio = (y, ynext) -> HMqTo(ynext, ynext) - HMqTo(y, y)
    Hia = (y, ynext) -> HMqTa(ynext, ynext) - HMqTa(y, y)

    function Hstiffo(y)
        σo = y[1] * scaleTraj
        return [0   σo   0   0   0]
    end
    function Hstiffa(y)
        σo = y[1] * scaleTraj
        return [0   0   0   m.ka*σo   m.ka*(-σo^3/3)]
    end
    function Hdamp(y)
        σ̇o = y[2] * scaleTraj * dtold
        return [0   0   σ̇o   0    0]
    end
    function Hvel(y)
        # Such that Hvel*pt[middle i.e./dt] = act. frame vel.
        φ = y[1]
        dφ = y[2] * dtold
        return [0   0   0    dφ   dφ*(-φ^2)]
    end
    
    if debugComponents
        return yo, Hio, Hia, Hstiffo, Hstiffa, Hdamp
    end
    # ----------------
    
    # Functions to output
    "Takes in a Δy in output coords"
    function Hk(k, Δyk, Δykp1)
        y = yo(k) + Δyk
        ynext = yo(k+1) + Δykp1
        # # This is inefficient since dydt is being called twice but fix later TODO:
        # yc, dyc = cu.collocationStates(m, opt, y, ynext, uk(k), uk(min(k+1,N)), param, dtold)
        H_dt2 = HMqT(y, ynext) - HMqT(y, y)
        H_dt1 = Hdamp(y)
        H_dt0 = Hstiffo(y) + Hstiffa(y)
        # With new nonlinear transmission need to break apart H
        σo = y[1]
        return [cu.Hτ(H_dt2, σo)  cu.Hτ(H_dt1, σo)   cu.Hτ(H_dt0, σo)], cu.Hτ(Hvel(y), σo)
    end
    
    # For a traj, H(yk, ykp1, Fk) * pb = B uk for each k
    B = reshape([1.0], 1, 1)

    return Hk, yo, uk, B, N
end

