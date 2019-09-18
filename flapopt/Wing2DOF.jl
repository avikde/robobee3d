
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
- R = wing length. This is irrelevant for the 2DOF model, but choosing one for now as 20mm. This is used to (a) calculate the aero force, and (b) in the plotting to show the "stroke angle" as a more intuitive version of σ.

NOTE:
- The "T" above is unitless. You can intuitively think of a wing "stroke angle" ~= T q[1] / (R / 2) (using σ as the arc length and "R/2" as the radius). This is distinct from the Toutput (in rad/m), and they are related as T ~= Toutput ⋅ (R / 2).
- Reasonable values: Toutput = 2666 rad/m in the current two-wing vehicle; 2150 rad/m in X-Wing; 3333 rad/m in Kevin Ma's older vehicles. With the R above, this suggests T ~= 20-30.
"""
struct Wing2DOFModel <: controlutils.Model
    # params
    mspar::Float64 # [mg]
    mwing::Float64 # [mg]
    kσ::Float64 # [mN/mm]
    bσ::Float64 # [mN/(mm/ms)]
    kΨ::Float64 # [mN-mm/rad]
    bΨ::Float64 # [mN-mm/(rad/ms)]
end

# Fixed params -----------------
const R = 20.0
const γ = 5.0 # wing shape fitting

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


"Returns paero [mm], Jaero, Faero [mN]"
function w2daero(y::AbstractArray, u::AbstractArray, _params::Vector)
    CLmax = 1.8
    CDmax = 3.4
    CD0 = 0.4
    ρ = 1.225e-3 # [mg/(mm^3)]
    
    # unpack
    cbar, T = _params
    σ, Ψ, σ̇, Ψ̇ = (@SVector [T, 1, T, 1]) .* y # [mm, rad, mm/ms, rad/ms]
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
    Faero = 1/2 * ρ * cbar * R * σ̇^2 * Caero * sign(-σ̇) # [mN]

    return paero, Jaero, Faero
end

"Continuous dynamics second order model"
function cu.dydt(m::Wing2DOFModel, y::AbstractArray, u::AbstractArray, _params::Vector)::AbstractArray
    # unpack
    cbar, T = _params
    # cbar, T, kσ = _params
    σ, Ψ, σ̇, Ψ̇ = (@SVector [T, 1, T, 1]) .* y # [mm, rad, mm/ms, rad/ms]
    # NOTE: for optimizing transmission ratio
    # Thinking of y = (sigma_actuator, psi, dsigma_actuator, dpsi)
    # u = (tau_actuator)
    # sigma = sigma_actuator * T; tau = tau_actuator / T
    cΨ = cos(Ψ)
    sΨ = sin(Ψ)

    Iwing = m.mwing * cbar^2 # cbar is in mm

    # inertial terms
    M = @SMatrix [m.mspar+m.mwing   cbar*m.mwing*cΨ; cbar*m.mwing*cΨ   Iwing+cbar^2*m.mwing]
    corgrav = @SVector [m.kσ*σ - cbar*m.mwing*sΨ*Ψ̇^2, m.kΨ*Ψ]
    # non-lagrangian terms
    τdamp = @SVector [-m.bσ * σ̇, -m.bΨ * Ψ̇]
    _, Jaero, Faero = w2daero(y, u, _params)
    τaero = Jaero' * Faero # units of [mN, mN-mm]
    # input
    τinp = @SVector [u[1]/T, 0]

    ddq = inv(M) * (-corgrav + τdamp + τaero + τinp)
    # return ddq
    return [y[3], y[4], ddq[1] / T, ddq[2]]
end

# Create an initial traj --------------

"Simulate with an OL signal at a freq"
function openloopResponse(m::Wing2DOFModel, opt::cu.OptOptions, freq::Real, params::Vector)
    # Create a traj
    σmax = cu.limits(m)[end][1]
    tend = 100.0 # [ms]
    function controller(y, t)
        return 50.0 * sin(2*π*freq*t)
    end
    vf(y, p, t) = cu.dydt(m, y, [controller(y, t)], params)
    # OL traj1
    simdt = 0.1 # [ms]
    teval = collect(0:simdt:tend) # [ms]
    prob = ODEProblem(vf, [0.2,0.,0.,0.], (teval[1], teval[end]))
    sol = solve(prob, saveat=teval)
    
    # Deduce metrics
    t = sol.t[end-100:end]
    y = hcat(sol.u[end-100:end]...)
    σmag = (maximum(y[1,:]) - minimum(y[1,:])) * params[2] / (R/2)
    Ψmag = maximum(y[2,:]) - minimum(y[2,:])
    # relative phase? Hilbert transform?

    # # Plot
    # σt = plot(sol, vars=3, ylabel="act vel [m/s]")
    # Ψt = plot(sol, vars=2, ylabel="hinge ang [r]")
    # plot(σt, Ψt, layout=(2,1))
    # gui()

    return [σmag, Ψmag]
end

"""
freq [kHz]; posGains [mN/mm, mN/(mm-ms)]; [mm, 1]
Example: trajt, traj0 = Wing2DOF.createInitialTraj(0.15, [1e3, 1e2], params0)
"""
function createInitialTraj(m::Wing2DOFModel, opt::cu.OptOptions, N::Int, freq::Real, posGains::Vector, params::Vector)
    # Create a traj
    σmax = cu.limits(m)[end][1]
    tend = 100.0 # [ms]
    function controller(y, t)
        # Stroke pos control
        σdes = 0.9 * σmax * sin(freq * 2 * π * t)
        return posGains[1] * (σdes - y[1]) - posGains[2] * y[3]
    end
    vf(y, p, t) = cu.dydt(m, y, [controller(y, t)], params)
    # OL traj1
    simdt = 0.1 # [ms]
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

    starti = 170
    olRange = starti:2:(starti + 2*N)
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
    cbar, T = params[1]
    # If plot is given a matrix each column becomes a different line
    σt = plot(t, hcat([traj[@view liy[1,:]] * T / (R/2) for traj in trajs]...), marker=:auto, ylabel="stroke ang [r]", title="δt=$(round(δt; sigdigits=4))ms; c=$(round(cbar; sigdigits=4))mm, T=$(round(T; sigdigits=4))")
    
    Ψt = plot(t, hcat([traj[@view liy[2,:]] for traj in trajs]...), marker=:auto, legend=false, ylabel="hinge ang [r]")
    
    ut = plot(t, hcat([[traj[@view liu[1,:]];NaN] for traj in trajs]...), marker=:auto, legend=false, ylabel="stroke force [mN]")
    Nt = length(trajs)
    # Plot of aero forces at each instant
    function aeroPlotVec(_traj::Vector, _param)
        Faerok = k -> w2daero(_traj[@view liy[:,k]], _traj[@view liu[:,k]], _param)[end]
        Faeros = hcat([Faerok(k) for k=1:N]...)
        return [Faeros[2,:]' NaN]'
    end
    aerot = plot(t, hcat([aeroPlotVec(trajs[i], params[i]) for i=1:Nt]...), marker=:auto, legend=false, ylabel="lift [mN]")
    # Combine the subplots
	return (σt, Ψt, ut, aerot)
end

function drawFrame(m::Wing2DOFModel, yk, uk, param; Faeroscale=1.0)
    cbar, T = param
    paero, _, Faero = w2daero(yk, uk, param)
    wing1 = [yk[1] * T;0] # wing tip
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

function plotParams(m::Wing2DOFModel, opt::cu.OptOptions, traj::Vector, args...; μ::Float64=1e-1)
    # First plot the param landscape
    p1 = 0:0.5:5.0
    p2 = 5:2:40

    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    Ng = opt.boundaryConstraint == cu.SYMMETRIC ? (N+2)*ny : (N+1)*ny
    g = zeros(Ng)#cu.gbounds(m, opt, traj)[1]
    f(p1, p2) = begin
        cu.gvalues!(g, m, opt, traj, [p1,p2], traj[1:4])
        cu.Jobj(m, opt, traj, [p1,p2]) + μ/2 * g' * g
    end

    paramLandscape = contour(p1, p2, f, fill=true)

    # Now plot the path taken
    params = hcat(args...) # Np x Nsteps
    plot!(paramLandscape, params[1,:], params[2,:], marker=:auto)
    
    return (paramLandscape)
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
function cu.robj(m::Wing2DOFModel, opt::cu.OptOptions, traj::AbstractArray, params::AbstractArray)::AbstractArray
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    
	yk = k -> @view traj[liy[:,k]]
	uk = k -> @view traj[liu[:,k]]
    
    Favg = @SVector zeros(2)
    for k = 1:N
        paero, _, Faero = w2daero(yk(k), uk(k), params)
        Favg += Faero
    end
    # Divide by the total time
    Favg /= (N) # [mN]
    # avg lift
    return [Favg[2] - 100]
end

# param opt stuff ------------------------

function cu.paramAffine(m::Wing2DOFModel, opt::cu.OptOptions, traj::AbstractArray, param::AbstractArray, Ruu::AbstractArray, Ryu::AbstractArray)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)

    # Make a new traj where the dynamics constraint is satisfied exactly
    traj1 = copy(traj)
	yk = k -> @view traj1[liy[:,k]]
    uk = k -> @view traj1[liu[:,k]]
    for k=1:N
        traj1[liy[:,k+1]] = yk(k) + δt * cu.dydt(m, yk(k), uk(k), param)
    end

    # Param stuff
    cbar, T = param
    # lumped parameter vector
    # pb = [T, T^2, cbar*T, cbar*T^2, cbar^2*T] # NOTE the actual param values are only needed for the test mode
    pb = [T^2, cbar*T]
    npb = length(pb)
    # This multiplies pbar from the left to produce the right side
    Iwing = m.mwing * cbar^2 # cbar is in mm
    
    # If test is true, it will test the affine relation
    test = false

    function HMq(ypos, yvel)
        σa, Ψ, σ̇adum, Ψ̇dum = ypos
        σadum, Ψdum, σ̇a, Ψ̇ = yvel
        return [0   σ̇a*(m.mspar+m.mwing)   Ψ̇*m.mwing*cos(Ψ)   0   0;
        Ψ̇*Iwing   0   0   σ̇a*m.mwing*cos(Ψ)   Ψ̇*m.mwing]
    end

    function HCgJ(y, F)
        σa, Ψ, σ̇a, Ψ̇ = y
        # See notes: this F stuff is w2d specific
        Ftil = -F/cbar
        rcop = 0.25 + 0.25 / (1 + exp(5.0*(1.0 - 4*(π/2 - abs(Ψ))/π))) # [(6), Chen (IROS2016)]

        # FIXME: this orig version probably has a negative sign error on the dynamics terms
        # return [0   -m.kσ*σa-m.bσ*σ̇a   Ftil[1]   0   0;
        # -m.kΨ*Ψ-m.bΨ*Ψ̇   0   0   -σ̇a*Ψ̇*m.mwing*sin(Ψ)   rcopnondim*(Ftil[1]*cos(Ψ) + Ftil[2]*sin(Ψ))]
        # 1st order version TODO: check the J^T*F
        return [0   m.kσ*σa+m.bσ*σ̇a   -Ψ̇^2*m.mwing*sin(Ψ)+Ftil[1]   0   0;
        m.kΨ*Ψ+m.bΨ*Ψ̇   0   0   0   rcop*(Ftil[1]*cos(Ψ) + Ftil[2]*sin(Ψ))]
    end

    # this is OK - see notes
    # Htil = (y, ynext, F) -> HMq(ynext) - HMq(y) + δt*HCgJ(y, F)
    # 1st order integration mods
    Htil = (y, ynext, F) -> HMq(y, ynext) - HMq(y, y) + δt*HCgJ(y, F)

    # For a traj, H(yk, ykp1, Fk) * pb = B uk for each k
    B = [1.0, 0.0]
    Bdag = (B' * B) \ B'
    P = zeros(npb, npb)
    q = zeros(npb)
    if test
        Hpb = zeros(ny÷2, N)
        Bu = similar(Hpb)
        errk = zeros(ny, N)
    end

    for k=1:N
        _, Jaero, Faero = w2daero(yk(k), uk(k), param)
        # Add same F as the traj produced NOTE: this has the assumption that the interaction force is held constant.
        Hh = Htil(yk(k), yk(k+1), Faero)[:,2:3] # FIXME: other cols are zero
        P += Hh' * Bdag' * Ruu * Bdag * Hh
        q += Hh' * Bdag' * Ryu' * yk(k)
        if test
            errk[:,k] = -yk(k+1) + yk(k) + δt * cu.dydt(m, yk(k), uk(k), param)
            Hpb[:,k] = Hh * pb
            Bu[:,k] = δt * B * uk(k)[1]
        end
    end
    if test
        return Hpb, Bu, errk
    else
        return Symmetric(P), q
    end
end


