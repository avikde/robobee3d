
using LinearAlgebra, StaticArrays, DifferentialEquations
using Plots; gr()

import controlutils
cu = controlutils

"""
2DOF wing model config:
- q[1] = actuator displacement (in mm)
- q[2] = hinge angle (in rad)
- y = (q, dq)
- u = [τ], where τ = actuator force in (mN)

Units:
- Instead of using SI units, I am using mm for length unit, mN for force unit, mg for mass, ms for time.
- with this, 1mN = 1(mg)*1(mm)/1(ms)^2 works out
- density 1kg/m^3 = 1e-3 mg/mm^3
- rotational inertia 1kg*m^2 = 1e9 mg*mm^2
- 1 N-m/rad = 1e6 mN-mm/rad
- 1 N-m/(rad/s) = 1e9 mN-mm/(rad/ms)

Reasonable limits:
- MRL Bee: umax = 75mN, q[1]max = 300 μm, mass = 25 mg
- MRL HAMR: umax = 270mN, q[1]max = 247 μm, mass = 100 mg [Doshi et. al. (2015)]

Params:
- cbar = wing chord
- T = transmission ratio. i.e. wing spar displacement σ = T q[1].
- R = wing length. This is irrelevant for the 2DOF model, but choosing one for now from [Jafferis (2016)]. This is used to (a) calculate the aero force, and (b) in the plotting to show the "stroke angle" as a more intuitive version of σ.

NOTE:
- The "T" above is unitless. You can intuitively think of a wing "stroke angle" ~= T q[1] / (R / 2) (using σ as the arc length and "R/2" as the radius). This is distinct from the Toutput (in rad/m), and they are related as T ~= Toutput ⋅ (R / 2).
- Reasonable values: Toutput = 2666 rad/m in the current two-wing vehicle; 2150 rad/m in X-Wing; 3333 rad/m in Kevin Ma's older vehicles. With the R above, this suggests T ~= 20-30.
"""
struct Wing2DOFModel <: controlutils.Model
    # overall scaling
    R::Float64 # [mm]
    # params
    kσ::Float64 # [mN/mm]
    bσ::Float64 # [mN/(mm/ms)]
    # Actuator params
    ma::Float64 # [mg]; effective
    ba::Float64 # [mN/(mm/ms)]
    ka::Float64 # [mN/mm]
    bCoriolis::Bool # true to include hinge->stroke Coriolis term or leave out for testing
end

# Fixed params -----------------
const γ = 0.5 # location of mwing lumped mass is γ*cbar down from the spar

function cu.dims(m::Wing2DOFModel)::Tuple{Int, Int}
    return 4, 1
end

function cu.pdims(m::Wing2DOFModel)::Int
    return 2
end

function cu.plimits(m::Wing2DOFModel)
    return [0.1, 1.0], [5.0, 100.0]
end

function cu.limits(m::Wing2DOFModel)::Tuple{Vector, Vector, Vector, Vector}
    # This is based on observing the OL trajectory. See note on units above.
    umax = @SVector [75.0] # [mN]
    umin = -umax
    xmax = @SVector [300e-3, 1.5, 100, 100] # [mm, rad, mm/ms, rad/ms]
    xmin = -xmax
    return umin, umax, xmin, xmax
end

function cu.limitsTimestep(m::Wing2DOFModel)::Tuple{Float64, Float64}
	return 0.01, 0.07
end


"Returns paero [mm], Jaero, Faero [mN]. Takes in y in *output coordinates*"
function w2daero(m::Wing2DOFModel, yo::AbstractArray, _params::Vector)
    CLmax = 1.8
    CDmax = 3.4
    CD0 = 0.4
    ρ = 1.225e-3 # [mg/(mm^3)]
    
    # unpack
    cbar, τ1, mwing, kΨ, bΨ, τ2 = _params
    σ, Ψ, σ̇, Ψ̇ = yo # [mm, rad, mm/ms, rad/ms]
    cΨ = cos(Ψ)
    sΨ = sin(Ψ)

    # CoP kinematics
    wing1 = @SVector [σ, 0]
    RΨ = @SMatrix [cΨ -sΨ; sΨ cΨ] # Ψ > 0 => hinge looks like \; Ψ < 0 => hinge looks like /
    α = π/2 - Ψ # AoA
    rcopnondim = 0.25 + 0.25 / (1 + exp(5.0*(1.0 - 4*(π/2 - abs(Ψ))/π))) # [(6), Chen (IROS2016)]
    paero = wing1 + RΨ * @SVector [0, -rcopnondim*cbar]
    # approx don't include the variation in COP in this
    Jaero = @SMatrix [1 rcopnondim*cbar * cΨ; 0 rcopnondim*cbar * sΨ]
    
    # Aero force
    #=
    Use this in Mathematica to debug the signs.
    Manipulate[
    \[Alpha] = -\[Pi]/2 + \[Psi];
    Graphics@
    Arrow[{{0, 0},
        {((CDmax + CD0)/2 - (CDmax - CD0)/2*Cos[2 \[Alpha]]), CLmax Sin[2 \[Alpha]]} Sign[-d\[Sigma]]}],
    {d\[Sigma], -1, 1}, {\[Psi], -\[Pi]/2, \[Pi]/2}]
    =#
    Caero = @SVector [((CDmax + CD0)/2 - (CDmax - CD0)/2 * cos(2α)), CLmax * sin(2α)]
    Faero = 1/2 * ρ * cbar * m.R * σ̇^2 * Caero * sign(-σ̇) # [mN]

    return paero, Jaero, Faero
end

"Continuous dynamics second order model"
function cu.dydt(m::Wing2DOFModel, y::AbstractArray, u::AbstractArray, _params::Vector)::AbstractArray
    # unpack
    cbar, τ1, mwing, kΨ, bΨ, τ2 = _params
    yo, T, τfun, τifun = cu.transmission(m, y, _params)
    σ, Ψ, σ̇, Ψ̇ = yo
    # NOTE: for optimizing transmission ratio
    # Thinking of y = (sigma_actuator, psi, dsigma_actuator, dpsi)
    # u = (tau_actuator)
    # sigma = sigma_actuator * T; tau = tau_actuator / T
    cΨ = cos(Ψ)
    sΨ = sin(Ψ)

    # inertial terms
    M = @SMatrix [mwing + m.ma/T^2   γ*cbar*mwing*cΨ;  γ*cbar*mwing*cΨ   2*cbar^2*γ^2*mwing]
    corgrav = @SVector [(m.kσ*σ + m.ka/T*τifun(σ)) + (m.bCoriolis ? -γ*cbar*mwing*sΨ*Ψ̇^2 : 0), kΨ*Ψ]
    # non-lagrangian terms
    τdamp = @SVector [-(m.bσ + m.ba/T^2) * σ̇, -bΨ * Ψ̇]
    _, Jaero, Faero = w2daero(m, yo, _params)
    τaero = Jaero' * Faero # units of [mN, mN-mm]
    # input
    τinp = @SVector [u[1]/T, 0]

    ddq = inv(M) * (-corgrav + τdamp + τaero + τinp)
    # return ddq
    return [y[3], y[4], ddq[1] / T, ddq[2]]
end

# Create an initial traj --------------

# "Simulate with an OL signal at a freq"
# function openloopResponse(m::Wing2DOFModel, opt::cu.OptOptions, freq::Real, params::Vector)
#     # Create a traj
#     σmax = cu.limits(m)[end][1]
#     tend = 100.0 # [ms]
#     function controller(y, t)
#         return 50.0 * sin(2*π*freq*t)
#     end
#     vf(y, p, t) = cu.dydt(m, y, [controller(y, t)], params)
#     # OL traj1
#     simdt = 0.1 # [ms]
#     teval = collect(0:simdt:tend) # [ms]
#     prob = ODEProblem(vf, [0.2,0.,0.,0.], (teval[1], teval[end]))
#     sol = solve(prob, saveat=teval)
    
#     # Deduce metrics
#     t = sol.t[end-100:end]
#     y = hcat(sol.u[end-100:end]...)
#     σmag = (maximum(y[1,:]) - minimum(y[1,:])) * params[2] / (m.R/2)
#     Ψmag = maximum(y[2,:]) - minimum(y[2,:])
#     # relative phase? Hilbert transform?

#     # # Plot
#     # σt = plot(sol, vars=3, ylabel="act vel [m/s]")
#     # Ψt = plot(sol, vars=2, ylabel="hinge ang [r]")
#     # plot(σt, Ψt, layout=(2,1))
#     # gui()

#     return [σmag, Ψmag]
# end

"""
freq [kHz]; posGains [mN/mm, mN/(mm-ms)]; [mm, 1]
Example: trajt, traj0 = Wing2DOF.createInitialTraj(0.15, [1e3, 1e2], params0)
"""
function createInitialTraj(m::Wing2DOFModel, opt::cu.OptOptions, N::Int, freq::Real, posGains::Vector, params::Vector, starti)
    # Create a traj
    σmax = 0.2#cu.limits(m)[end][1]
    tend = 100.0 # [ms]
    function controller(y, t)
        # Stroke pos control
        σdes = σmax * sin(freq * 2 * π * t)
        return posGains[1] * (σdes - y[1]) - posGains[2] * y[3]
    end
    vf(y, p, t) = cu.dydt(m, y, [controller(y, t)], params)
    # OL traj1
    simdt = opt.fixedδt # [ms]
    teval = collect(0:simdt:tend) # [ms]
    prob = ODEProblem(vf, [0.2,0.,0.,0.], (teval[1], teval[end]))
    sol = solve(prob, saveat=teval)

    # # Animate whole traj
    # Nt = length(sol.t)
    # @gif for k = 1:3:Nt
    #     yk = sol.u[k]
    #     uk = [controller(yk, sol.t[k])]
    #     drawFrame(m, yk, uk, params)
    # end
    # # Plot
    # σt = plot(sol, vars=3, ylabel="act vel [m/s]")
    # Ψt = plot(sol, vars=2, ylabel="hinge ang [r]")
    # plot(σt, Ψt, layout=(2,1))
    # gui()

    # expectedInterval = opt.boundaryConstraint == cu.SYMMETRIC ? 1/(2*freq) : 1/freq # [ms]
    # expectedPts = expectedInterval / simdt

    olRange = starti:(starti + N)
    trajt = sol.t[olRange]
    δt = trajt[2] - trajt[1]
    olTrajaa = sol.u[olRange] # 23-element Array{Array{Float64,1},1} (array of arrays)
    olTraju = [controller(olTrajaa[i], trajt[i]) for i in 1:N] # get u1,...,uN
    traj0 = [vcat(olTrajaa...); olTraju] # dirtran form {x1,..,x(N+1),u1,...,u(N),δt}
    if opt.vart
        traj0 = [traj0; δt]
    else
        println("Initial traj δt=", δt, ", opt.fixedδt=", opt.fixedδt)
    end

    return trajt .- trajt[1], traj0
end

function plotTrajs(m::Wing2DOFModel, opt::cu.OptOptions, t::Vector, params, trajs)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, trajs[1])
	Ny = (N+1)*ny
    # stroke "angle" = T*y[1] / R
    function strokeAng(traj, param)
        σo = [cu.transmission(m, traj[@view liy[:,k]], param)[1][1] for k=1:N+1]
        return σo / (m.R/2)
    end
    # If plot is given a matrix each column becomes a different line
    σt = plot(t, hcat([strokeAng(trajs[i], params[i]) for i in 1:length(trajs)]...), linewidth=2, ylabel="stroke ang [r]")# title="δt=$(round(δt; sigdigits=4))ms; c=$(round(cbar; sigdigits=4))mm, T=$(round(T; sigdigits=4))")
    
    Ψt = plot(t, hcat([traj[@view liy[2,:]] for traj in trajs]...), linewidth=2, legend=false, ylabel="hinge ang [r]")
    
    ut = plot(t, hcat([[traj[@view liu[1,:]];NaN] for traj in trajs]...), linewidth=2, legend=false, ylabel="stroke force [mN]")
    Nt = length(trajs)
    # Plot of aero forces at each instant
    function aeroPlotVec(_traj::Vector, _param, ind)
        cbar, τ1, mwing, kΨ, bΨ, τ2 = _param
        Faerok = k -> w2daero(m, cu.transmission(m, _traj[@view liy[:,k]], _param)[1], _param)[end] * 1000 / 9.81 # to mg
        Faeros = hcat([Faerok(k) for k=1:N]...)
        return [Faeros[ind,:]' NaN]'
    end
    liftt = plot(t, hcat([aeroPlotVec(trajs[i], params[i], 2) for i=1:Nt]...), linewidth=2, legend=false, ylabel="lift [mg]")
    dragt = plot(t, hcat([aeroPlotVec(trajs[i], params[i], 1) for i=1:Nt]...), linewidth=2, legend=false, ylabel="drag [mg]")
    # Combine the subplots
	return (σt, Ψt, ut, liftt, dragt)
end

function plotParamImprovement(m::Wing2DOFModel, opt::cu.OptOptions, t::Vector, params, trajs, paramObj::Function)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, trajs[1])

    # The param space plots
    pls = plotParams(m, opt, trajs[1], paramObj, params...)
    # Traj plots
    σt, Ψt, ut, liftt, dragt = plotTrajs(m, opt, t, params, trajs)

    return pls..., σt, Ψt, ut, liftt, dragt
end

function compareTrajToDAQ(m::Wing2DOFModel, opt::cu.OptOptions, t::Vector, param, traj, lift, drag)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
	Ny = (N+1)*ny
	# stroke "angle" = T*y[1] / R
    cbar, τ1, mwing, kΨ, bΨ, τ2, dt = param
    # Plot of aero forces at each instant
    function aeroPlotVec(_traj::Vector, _param, i)
        Faerok = k -> w2daero(m, _traj[@view liy[:,k]], _param)[end]
        Faeros = hcat([Faerok(k) for k=1:N]...)
        return [Faeros[i,:]' NaN]'
    end
    liftp = plot(t, lift, marker=:auto, legend=false, label="daq", ylabel="lift [mN]")
    plot!(liftp, t, aeroPlotVec(traj, param, 2), marker=:auto, legend=false, label="pred")
    dragp = plot(t, drag, marker=:auto, label="daq", ylabel="drag [mN]")
    plot!(dragp, t, aeroPlotVec(traj, param, 1), marker=:auto, label="pred")

    # Get the basic kinematics
    σt, Ψt, ut, _ = plotTrajs(m, opt, t, [param], [traj])

    # Combine the subplots
	return (σt, Ψt, ut, liftp, dragp)
end

function drawFrame(m::Wing2DOFModel, yk, uk, param; Faeroscale=1.0)
    cbar, τ1, mwing, kΨ, bΨ, τ2, dt = param
    yo, T, _, _ = cu.transmission(m, yk, param)
    paero, _, Faero = w2daero(m, yo, param)
    wing1 = [yo[1];0] # wing tip
    wing2 = wing1 + normalize(paero - wing1)*cbar
    # draw wing
    w = plot([wing1[1]; wing2[1]], [wing1[2]; wing2[2]], color=:gray, aspect_ratio=:equal, linewidth=5, legend=false, xlims=(-15,15), ylims=(-3,3))
    # Faero
    FaeroEnd = paero + Faeroscale * Faero
    plot!(w, [paero[1], FaeroEnd[1]], [paero[2], FaeroEnd[2]], color=:red, linewidth=2, line=:arrow)
    # stroke plane
    plot!(w, [-T,T], [0, 0], color=:black, linestyle=:dash)
    return w
end

function animateTrajs(m::Wing2DOFModel, opt::cu.OptOptions, params, trajs)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, trajs[1])

    # TODO: simulate with a waypoint controller for more steps
    
	yk(traj, k) = @view traj[liy[:,k]]
    uk(traj, k) = @view traj[liu[:,k]]
    
    Nanim = opt.boundaryConstraint == :symmetric ? 2*N : N
    Nt = length(trajs)

    function drawTrajFrame(traj, param, k)
        _yk = opt.boundaryConstraint == :symmetric && k > N ? -yk(traj, k-N) : yk(traj, k)
        _uk = opt.boundaryConstraint == :symmetric && k > N ? -uk(traj, k-N) : uk(traj, k)
        return drawFrame(m, _yk, _uk, param)
    end
    @gif for k=1:Nanim
        plot([drawTrajFrame(trajs[i], params[i], k) for i in 1:Nt]..., layout=(Nt,1))
    end

    # return wingdraw 
end

function plotParams(m::Wing2DOFModel, opt::cu.OptOptions, traj::Vector, paramObj::Function, args...; μ::Float64=1e-1, Vclip=50000)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    # First plot the param landscape
    pranges = [
        0:0.25:6.0, # cbars
        10.0:1.0:50, # Ts
        0.1:0.1:3.0, # mwings
        0.1:1.0:20.0, # kΨs
        0.1:1.0:20.0 # bΨs
    ]
    labels = [
        "chord",
        "T",
        "mwing",
        "hinge k",
        "hinge b"
    ]

    # different param vectors passed in
    params = hcat(args...) # Np x Nsteps
    param0 = args[end] # for the slices use the last param

    # # Old: f defined here
    # Ng = opt.boundaryConstraint == cu.SYMMETRIC ? (N+2)*ny : (N+1)*ny
    # g = zeros(Ng)#cu.gbounds(m, opt, traj)[1]
    # f(p1, p2) = begin
    #     cu.gvalues!(g, m, opt, traj, [p1,p2], traj[1:4])
    #     cu.Jobj(m, opt, traj, [p1,p2]) + μ/2 * g' * g
    # end

    function plotSlice(i1, i2)
        function f(p1, p2)
            parg = copy(param0)
            parg[i1] = p1
            parg[i2] = p2
            V = paramObj(parg)
            return V > Vclip ? NaN : V
        end
        pl = contour(pranges[i1], pranges[i2], f, fill=true, seriescolor=cgrad(:bluesreds), xlabel=labels[i1], ylabel=labels[i2])
        # Now plot the path taken
        plot!(pl, params[i1,:], params[i2,:], marker=:auto, legend=false)
        # just in case
        xlims!(pl, (pranges[i1][1], pranges[i1][end]))
        ylims!(pl, (pranges[i2][1], pranges[i2][end]))
        return pl
    end
    
    return (plotSlice(1, 2), # cbar, T slice
        plotSlice(1, 3), # cbar, mwing slice
        plotSlice(2, 3), # T, mwing slice
        plotSlice(4, 5) # T, khinge slice
    )
end

# Test applying euler integration to the initial traj
function eulerIntegrate(m::cu.Model, opt::cu.OptOptions, traj::Vector, params::Vector)
    trajei = copy(traj)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)

    yk(k) = @view trajei[liy[:,k]]
    uk(k) = @view trajei[liu[:,k]]
    for k=1:N
        trajei[liy[:,k+1]] = cu.ddynamics(m, yk(k), uk(k), params, δt)
    end
    return trajei
end

# "Cost function components" ------------------

"Objective to minimize"
function cu.robj(m::Wing2DOFModel, opt::cu.OptOptions, traj::AbstractArray, param::AbstractArray)::AbstractArray
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    
	yk = k -> @view traj[liy[:,k]]
	uk = k -> @view traj[liu[:,k]]
    
    cbar, τ1, mwing, kΨ, bΨ, τ2, dt = param

    Favg = @SVector zeros(2)
    for k = 1:N
        paero, _, Faero = w2daero(m, cu.transmission(m, yk(k), param)[1], param)
        Favg += Faero
    end
    # Divide by the total time
    Favg /= (N) # [mN]
    # avg lift
    return [Favg[2] - 100]
end

# param opt stuff ------------------------

function cu.transmission(m::Wing2DOFModel, y::AbstractArray, _param::Vector; o2a=false)
    cbar, τ1, mwing, kΨ, bΨ, τ2 = _param
    τfun = σa -> τ1*σa + τ2/3*σa^3
    # Series from Mathematica
    τifun = σo -> σo/τ1 - τ2*σo^3/(3*τ1^4) # + O[σo^4]
    if !o2a
        Dτfun = σa -> τ1 + τ2*σa^2
        T = Dτfun(y[1]) # "gear ratio"

        y2 = [τfun(y[1]), y[2], T*y[3], y[4]] # [mm, rad, mm/ms, rad/ms]
    else
        σo = y[1]
        T = τ1 + σo^2*τ2/τ1^2
        y2 = [τifun(σo), y[2], y[3]/T, y[4]] # [mm, rad, mm/ms, rad/ms]
    end
    return y2, T, τfun, τifun
end

function cu.paramLumped(m::Wing2DOFModel, param::AbstractArray)
    cbar, τ1, mwing, kΨ, bΨ, τ2, dt = param
    Tarr = [τ1, τ2]
    return [1, kΨ, bΨ, cbar, cbar^2, mwing, mwing*cbar, mwing*cbar^2], Tarr, dt
end

function cu.paramAffine(m::Wing2DOFModel, opt::cu.OptOptions, traj::AbstractArray, param::AbstractArray, POPTS::cu.ParamOptOpts, scaleTraj=1.0; debugComponents::Bool=false)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)

    yk = k -> @view traj[liy[:,k]]
    uk = k -> @view traj[liu[:,k]]

    # Param stuff
    cbar, τ1, mwing, kΨ, bΨ, τ2, dt = param
    # Need the original T to use output coords
    yo = k -> cu.transmission(m, yk(k), param)[1]

    # THESE FUNCTIONS USE OUTPUT COORDS -------------
    """Inertial"""
    function HMqTWithoutCoupling(ypos, yvel)
        σo, Ψ, σ̇odum, Ψ̇dum = ypos
        σodum, Ψdum, σ̇o, Ψ̇ = yvel
        return [0   0   0   0   0   σ̇o   0   0   σ̇o*m.ma   σ̇o*m.ma*(-σo^2);
        0   0   0   0   0   0   0   2*Ψ̇*γ^2    0   0]
    end
    function HMqTCoupling(ypos, yvel)
        σo, Ψ, σ̇odum, Ψ̇dum = ypos
        σodum, Ψdum, σ̇o, Ψ̇ = yvel
        return [0   0   0   0   0   0   γ*Ψ̇*cos(Ψ)   0   0   0;
        0   0   0   0   0   0   γ*σ̇o*cos(Ψ)   0    0   0]
    end
    HMqT(ypos, yvel) = HMqTWithoutCoupling(ypos, yvel) + HMqTCoupling(ypos, yvel)

    """Coriolis"""
    function HC(y)
        σo, Ψ, σ̇o, Ψ̇ = y
        return [0   0   0   0   0   0   (m.bCoriolis ? -γ*Ψ̇^2*sin(Ψ) : 0)   0    0   0;
        0   0   0    0   0   0   0   0   0   0]
    end
    """Stiffness/damping output"""
    function Hg(y)
        σo, Ψ, σ̇o, Ψ̇ = y
        return [m.kσ*σo+m.bσ*σ̇o   0   0   0   0   0   0   0    0   0;
        0   Ψ   Ψ̇    0   0   0   0   0   0   0]
    end
    """Stiffness/damping actuator"""
    function Hgact(y)
        σo, Ψ, σ̇o, Ψ̇ = y
        return [0   0   0   0   0   0   0   0    m.ba*σ̇o+m.ka*σo   m.ba*σ̇o*(-σo^2)+m.ka*(-σo^3/3);
        0   0   0    0   0   0   0   0   0   0]
    end
    """Ext force (aero)"""
    function HF(y, F)
        σo, Ψ, σ̇o, Ψ̇ = y
        # See notes: this F stuff is w2d specific
        rcop = 0.25 + 0.25 / (1 + exp(5.0*(1.0 - 4*(π/2 - abs(Ψ))/π))) # [(6), Chen (IROS2016)]

        if POPTS.Fext_pdep
            Ftil = -F/cbar

            # FIXME: this orig version probably has a negative sign error on the dynamics terms
            # return [0   -m.kσ*σa-m.bσ*σ̇a   Ftil[1]   0   0;
            # -m.kΨ*Ψ-m.bΨ*Ψ̇   0   0   -σ̇a*Ψ̇*m.mwing*sin(Ψ)   rcopnondim*(Ftil[1]*cos(Ψ) + Ftil[2]*sin(Ψ))]
            # 1st order version
            return [0   0   0   Ftil[1]   0   0   0   0   0   0;
            0   0   0    0   rcop*(Ftil[1]*cos(Ψ) + Ftil[2]*sin(Ψ))   0   0   0   0   0]
        else
            Ftil = -F
            # 1st order version
            return [Ftil[1]   0   0   0   0   0   0   0    0   0;
            0   0   0    rcop*(Ftil[1]*cos(Ψ) + Ftil[2]*sin(Ψ))   0   0   0   0   0   0]
        end
    end
    function HF2(y)
        _, _, Faero = w2daero(m, y, param)
        return HF(y, Faero)
    end

    if debugComponents
        return yo, HMqTWithoutCoupling, HMqTCoupling, HC, Hg, Hgact, HF2
    end

    HCgJT(y, F) = HC(y) + Hg(y) + Hgact(y) + HF(y, F)
    # ----------------

    # this is OK - see notes
    # Htil = (y, ynext, F) -> HMq(ynext) - HMq(y) + δt*HCgJ(y, F)
    # 1st order integration mods
    Htili = (y, ynext) -> HMqT(y, ynext) - HMqT(y, y)
    Htilni = HCgJT
    
    # Functions to output
    "Takes in a Δy in output coords"
    function Hk(k, Δyk, Δykp1)
        # Same Faero as before?
        _, _, Faero = w2daero(m, yo(k) + Δyk, param)
        # NOTE: it uses param *only for Faero*. Add same F as the traj produced
        Hi = Htili(yo(k) + Δyk, yo(k+1) + Δykp1)
        Hni = Htilni(yo(k) + Δyk, Faero)
        # Hh = Hi + δt*Hni # inertial and non-inertial components
        # With nonlinear transmission need to break apart H
        σo = (yo(k) + Δyk)[1]
        Hfortrans = Hh -> hcat(Hh[:,1:end-2], Hh[:,1:end-2]*σo^2, Hh[:,end-1:end])
        return [Hfortrans(Hi)  Hfortrans(Hni)]
    end
    
    # For a traj, H(yk, ykp1, Fk) * pb = B uk for each k
    B = reshape([1.0, 0.0], 2, 1)

    return Hk, yo, uk, B, N
end

function avgLift(m::Wing2DOFModel, opt::cu.OptOptions, traj::AbstractArray, param::AbstractArray)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    yk = k -> @view traj[liy[:,k]]
    uk = k -> @view traj[liu[:,k]]

    cbar, τ1, mwing, kΨ, bΨ, τ2, dt = param
    aa = 0
    for k=1:N
        _, _, Faero = w2daero(m, cu.transmission(m, yk(k), param)[1], param)
        aa += Faero[2] * δt
    end
    return aa / (N*δt)
end
