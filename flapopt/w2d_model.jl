
using Parameters, LinearAlgebra, DifferentialEquations, Dierckx, Statistics
using Plots

import controlutils
cu = controlutils

"""
2DOF wing model config:
- q[1] = stroke angle (in rad)
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
- T = transmission ratio
- Lw = wing length. This is irrelevant for the 2DOF model, but choosing one for now from [Jafferis (2016)]. This is used to (a) calculate the aero force, and (b) in the plotting to show the "stroke angle" as a more intuitive version of σ.

NOTE:
- Reasonable values: Toutput = 2666 rad/m in the current two-wing vehicle; 2150 rad/m in X-Wing; 3333 rad/m in Kevin Ma's older vehicles.
"""
@with_kw struct Wing2DOFModel <: controlutils.Model
    # params
    kbo::Array{Float64,1} = [0., 0.] # [mN/mm], [mN/(mm/ms)]
    # Actuator params
    ma::Float64 # [mg]; effective
    ba::Float64 = 0 # [mN/(mm/ms)]
    ka::Float64 # [mN/mm]
    bCoriolis::Bool = true # true to include hinge->stroke Coriolis term or leave out for testing
    r1h::Float64 = 0.49 # r1hat first moment -- see [Whitney (2010)]. Default from [Chen (2016)]
    r2h::Float64 = 0.551 # r2hat second moment -- see [Whitney (2010)]. Default from 0.929 * 0.49^0.732
    SEA::Bool = false # series-elastic https://github.com/avikde/robobee3d/pull/126
    kSEA::Float64 = 1000 # SEA spring const
    Amp::AbstractArray = [0.0,0.0] # Stroke, hinge amplitude for the desired output kinematics (set to 0 to not use)
end

# Fixed params -----------------
const γ = 0.5 # location of mwing lumped mass is γ*cbar down from the spar

function cu.dims(m::Wing2DOFModel)::Tuple{Int, Int}
    return m.SEA ? 6 : 4, 1
end

"Tried going directly from Doshi model-driven, but haven't been able to get that to match up"
function hingeParams(wΨ)
    # kΨ [mN-mm/rad], bΨ [mN-mm/(rad/ms)]
    return 3.0*wΨ, 1.8*wΨ
end

rcopnondim(Ψ) = 0.4#0.25 + 0.25 / (1 + exp(5.0*(1.0 - 4*(π/2 - abs(Ψ))/π))) # [(6), Chen (IROS2016)]

"Returns paero [mm], Jaero, Faero [mN]. Takes in y in *output coordinates*"
function w2daero(m::Wing2DOFModel, yo::AbstractArray, param::Vector)
    CLmax = 1.8
    CDmax = 3.4
    CD0 = 0.4
    ρ = 1.225e-3 # [mg/(mm^3)]
    
    # unpack
    cbar2, τ1, mwing, wΨ, τ2, Aw, dt = param
    AR = Aw/cbar2
    cbar = sqrt(cbar2)
    Lw = Aw/cbar
    ycp = Lw*m.r1h # approx as in [Chen (2016)]
    nq = length(yo)÷2
    φ, Ψ, dφ, dΨ = yo[[1,2,nq+1,nq+2]] # [rad, rad]

    # CoP kinematics
    # Ψ > 0 => hinge looks like \; Ψ < 0 => hinge looks like /
    α = π/2 - Ψ # AoA
    
    # From Mathematica
    rcnd = rcopnondim(Ψ)
    paero = [-ycp*sin(φ) - cbar*rcnd*cos(φ)*sin(Ψ), ycp*cos(φ) - cbar*rcnd*sin(φ)*sin(Ψ), -cbar*rcnd*cos(Ψ)]
    Jaero = [-ycp*cos(φ) + cbar*rcnd*sin(φ)*sin(Ψ)    -cbar*rcnd*cos(φ)*cos(Ψ);
            -ycp*sin(φ) - cbar*rcnd*cos(φ)*sin(Ψ)    -cbar*rcnd*cos(Ψ)*sin(φ);
            0     cbar*rcnd*sin(Ψ)]
    # drag/lift direction (normalized)
    eD = [dφ*cos(φ), dφ*sin(φ), 0]/(abs(dφ) + 1e-4) # so it is well-defined
    eL = [0, 0, 1]
    
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
    Caero = [((CDmax + CD0)/2 - (CDmax - CD0)/2 * cos(2α)), CLmax * sin(2α)]
    Faero = 1/2 * ρ * dφ^2 * Aw^2 * AR * m.r2h^2 * (Caero[1]*eD + Caero[2]*eL*sign(-dφ)) # [mN]

    return paero, Jaero, Faero
end

"Continuous dynamics second order model"
function cu.dydt(m::Wing2DOFModel, yo::AbstractArray, u::AbstractArray, param::Vector)::AbstractArray
    # unpack
    cbar2, τ1, mwing, wΨ, τ2, Aw, dt = param
    cbar = sqrt(cbar2)
    Lw = Aw/cbar
    ycp = Lw*m.r1h # approx as in [Chen (2016)]
    # These next two are "wingsubs" from the "Incorporate R" notes
    Izz = 0
    Ixx = cbar2*γ^2*mwing
    kΨ, bΨ = hingeParams(wΨ)
    nq = length(yo)÷2
    φ, Ψ, dφ, dΨ = yo[[1,2,nq+1,nq+2]] # [rad, rad]
    
    ya, T, τfun, τifun = cu.transmission(m, yo, param; o2a=true)

    # Lagrangian terms - from Mathematica
    M = [Izz + mwing*ycp^2 + 1/2*cbar2*mwing*γ^2*(1-cos(2*Ψ))   γ*cbar*mwing*ycp*cos(Ψ);  
        γ*cbar*mwing*ycp*cos(Ψ)     Ixx+cbar2*γ^2*mwing]
    cor1 = [cbar2*mwing*γ^2*sin(2*Ψ)*dφ*dΨ - γ*cbar*mwing*ycp*sin(Ψ)*dΨ^2, 
        -cbar2*mwing*γ^2*cos(Ψ)*sin(Ψ)*dφ^2]
    # NOTE: dropping τinv'' term
    corgrav = [m.kbo[1]*φ, kΨ*Ψ] + (m.bCoriolis ? cor1 : zeros(2))

    # non-lagrangian terms
    τdamp = [-(m.kbo[2] + m.ba/T^2) * dφ, -bΨ * dΨ]
    _, Jaero, Faero = w2daero(m, yo, param)
    τaero = Jaero' * Faero # units of [mN, mN-mm]
    # input
    τinp = [u[1]/T, 0]

    if m.SEA
        δa, dδa = yo[[3,nq+3]]
        # Add the new DOF for the actuator SEA-decoupled from the stroke
        M = vcat(hcat(M,zeros(2)),zeros(1,3))
        M[3,3] = m.ma
        seaForce = m.kSEA*(τfun(δa) - φ)
        append!(corgrav, T*seaForce + m.ka*δa)
        corgrav[1] -= seaForce
        append!(τdamp, 0)
        append!(τaero, 0)
        τinp = [0, 0, u[1]]
    else
        # Add reflected actuator properties
        M[1,1] += m.ma/T^2
        corgrav[1] += m.ka/T*τifun(φ)
    end

    ddq = inv(M) * (-corgrav + τdamp + τaero + τinp)
    # return ddq
    return [yo[nq+1:end]; ddq]
end

# Create an initial traj --------------

"""
freq [kHz]; posGains [mN/mm, mN/(mm-ms)]; [mm, 1]
Example: trajt, traj0 = Wing2DOF.createInitialTraj(0.15, [1e3, 1e2], params0)
"""
function createInitialTraj(m::Wing2DOFModel, opt::cu.OptOptions, N::Int, freq::Real, posGains::Vector, params::Vector, starti; uampl=65, h3=0.0, h2=0.0, posctrl=false, makeplot=false, trajstats=false, trajstats2=false, simdt=0.02, verbose=true, rawtraj=false, phaseoffs=0, dcoffs=0, powplot=false)
    # Create a traj
    φampl = 0.6 # output, only used if posctrl=true
    tend = 100.0 # [ms]
    function controller(y, t)
        if posctrl
            # Stroke pos control
            φdes = φampl * sin(freq * 2 * π * t + phaseoffs)
            dφdes = φampl * freq * 2 * π * cos(freq * 2 * π * t)
            dφ = m.SEA ? y[4] : y[3]
            return posGains[1] * (φdes - y[1]) + posGains[2] * (dφdes - dφ)
        else
            ph = freq * 2 * π * t + phaseoffs
            return uampl * ((1+h3)*sin(ph) + h3*sin(3*ph) + h2*sin(2*ph)) + dcoffs
        end
    end
    vf(y, p, t) = cu.dydt(m, y, [controller(y, t)], params)
    # OL traj1
    teval = collect(0:simdt:tend) # [ms]
    y0 = m.SEA ? [0.,0.01,0.,0.01,0.,0.] : [0.,0.01,0.01,0.]
    prob = ODEProblem(vf, y0, (teval[1], teval[end]))
    sol = solve(prob, saveat=teval)

    # # Animate whole traj
    # Nt = length(sol.t)
    # @gif for k = 1:3:Nt
    #     yk = sol.u[k]
    #     uk = [controller(yk, sol.t[k])]
    #     drawFrame(m, yk, uk, params)
    # end
    #
    function calcTrajAmplitudes()
        # don't return dirtran traj; only traj stats
        Nend = 1000
        solend = hcat(sol.u[end-Nend+1:end]...) # (ny,N) shaped
        vecrange(x) = maximum(x) - minimum(x)
        campl = k -> vecrange(solend[k,:])
        σa = [cu.transmission(m, solend[:,k], params; o2a=true)[1][1] for k=1:Nend]
        return [campl(1), campl(3), vecrange(σa)/2] # stroke, hinge ampl, act disp ampl
    end
    if trajstats
        return calcTrajAmplitudes()
    end
    if makeplot
        # Plot
        phit = plot(sol, vars=1, ylabel="stroke [rad]")
        dphit = plot(sol, vars=3, ylabel="stroke vel [rad/ms]")
        Ψt = plot(sol, vars=2, ylabel="hinge ang [r]")
        plot(phit, dphit, Ψt, layout=(3,1))
        gui()
        error("createInitialTraj")
    end
    if rawtraj
        return hcat(sol.u...)
    end
    # expectedInterval = opt.boundaryConstraint == cu.SYMMETRIC ? 1/(2*freq) : 1/freq # [ms]
    # expectedPts = expectedInterval / simdt

    # Get the traj timesteps to be correct using interp
    trajt = (starti:(starti + N))*opt.fixedδt
    soly = hcat(sol.u...) # ny,Nodesol
    spl = [Spline1D(sol.t, soly[i,:]) for i=1:ny]
    olTrajaa = hcat([spl[i](trajt) for i=1:ny]...) # ny arrays of size N each -> (N,ny)
    olTraju = [controller(olTrajaa[k,:], trajt[k]) for k in 1:N] # get u1,...,uN
    traj0 = [olTrajaa'[:]; olTraju] # dirtran form {x1,..,x(N+1),u1,...,u(N),δt}

    if trajstats2
        if powplot
	        ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj0)
            actvel = zeros(N)
            actu = zeros(N)
            for k=1:N
                T = cu.transmission(m, traj0[liy[:,k]], params; o2a=true)[2]
                actvel[k] = traj0[liy[3,k]]/T * actu[liu[1,k]]
            end
            pl = plot(trajt, actvel, label="actvel")
            plot!(pl, trajt, actu, label="actu")
            return pl
        end
        return [calcTrajAmplitudes(); avgLift(m, opt, traj0, params); trajMechPow(m, opt, traj0, params)]
    end

    # Some printing
    δt = trajt[2] - trajt[1]
    if opt.vart
        traj0 = [traj0; δt]
    elseif verbose
        println("Initial traj δt=", round(δt, digits=4), ", opt.fixedδt=", opt.fixedδt)
    end

    return trajt .- trajt[1], traj0
end

"Get actuator displacement from the output trajectory"
function actAng(m::Wing2DOFModel, opt::cu.OptOptions, traj, param)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)
    σa = [cu.transmission(m, traj[@view liy[:,k]], param; o2a=true)[1][1] for k=1:N+1]
    if m.SEA
        return [σa  traj[liy[3,:]]]
    else
        return σa
    end
end
"Get aero force components from the output trajectory. comp=:lift or :drag"
function trajAero(m::Wing2DOFModel, opt::cu.OptOptions, _traj::Vector, _param, comp; mg=true)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, _traj)
    cbar2, τ1, mwing, wΨ, τ2, Aw, dt  = _param
    Faerok = k -> w2daero(m, _traj[@view liy[:,k]], _param)[end] * (mg ? 1000 / 9.81 : 1) # to mg
    Faeros = hcat([Faerok(k) for k=1:N]...)
    trajAeroInd(ind) = [Faeros[ind,:];NaN] # need one more element
    
    return comp==:lift ? trajAeroInd(3) : sqrt.(trajAeroInd(1).^2 + trajAeroInd(2).^2) # drag has x,y components
end

function plotTrajs(m::Wing2DOFModel, opt::cu.OptOptions, params, trajs; legends=true)
	ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, trajs[1])
    Ny = (N+1)*ny
    nq = ny÷2
    Nt = length(trajs)

    # Empty versions of all the subplots
    stroket = plot(ylabel="stroke ang [deg]", ylims=(-60,60), legend=legends)# title="δt=$(round(δt; sigdigits=4))ms; c=$(round(cbar; sigdigits=4))mm, T=$(round(T; sigdigits=4))")
    Ψt = plot(ylabel="hinge ang [deg]", ylims=(-80,80), legend=legends)
    ut = plot(ylabel="stroke force [mN]", legend=legends)
    liftt = plot(ylabel="lift [mg]", legend=legends)
    dragt = plot(ylabel="drag [mg]", legend=legends)
    actt = plot(ylabel="act disp [mm]", #= ylims=(-0.3,0.3), =# legend=legends)
    phaset = plot(ylabel="Phase offs", legend=legends)
    for i=1:Nt
        traj, param = trajs[i], params[i]
        dt = param[end]
        t = 0:dt:(N)*dt
        plot!(stroket, t, rad2deg.(traj[@view liy[1,:]]), linewidth=2)
        plot!(actt, t, actAng(m, opt, traj, param), linewidth=2)
        plot!(Ψt, t, rad2deg.(traj[@view liy[2,:]]), linewidth=2)
        plot!(ut, t, [traj[@view liu[1,:]];NaN], linewidth=2)
        plot!(liftt, t, trajAero(m, opt, traj, param, :lift), linewidth=2)
        plot!(dragt, t, trajAero(m, opt, traj, param, :drag), linewidth=2)
        plot!(phaset, traj[@view liy[nq+1,:]], traj[@view liy[2,:]], linewidth=2)
    end
    # Return tuple
	return (stroket, Ψt, ut, liftt, dragt, actt, phaset)
end

function drawFrame(m::Wing2DOFModel, yk, uk, param; Faeroscale=1.0)
    cbar2, τ1, mwing, wΨ, τ2, Aw, dt = param
    yo, T, _, _ = cu.transmission(m, yk, param)
    paero, _, Faero = w2daero(m, yo, param)
    wing1 = [yo[1];0] # wing tip
    wing2 = wing1 + normalize(paero - wing1)*sqrt(cbar2)
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

# param opt stuff ------------------------

function cu.transmission(m::Wing2DOFModel, y::AbstractArray, _param::Vector; o2a=false)
    cbar2, τ1, mwing, wΨ, τ2, Aw, dt = _param
    τfun = σa -> τ1*σa + τ2/3*σa^3
    # Series from Mathematica
    τifun = φo -> φo/τ1 - τ2*φo^3/(3*τ1^4) # + O[φo^4]
    if !o2a
        Dτfun = σa -> τ1 + τ2*σa^2
        T = Dτfun(y[1]) # "gear ratio"

        y2 = [τfun(y[1]), y[2], T*y[3], y[4]] # [mm, rad, mm/ms, rad/ms]
    else
        φo = y[1]
        T = τ1 + φo^2*τ2/τ1^2
        y2 = [τifun(φo), y[2], y[3]/T, y[4]] # [mm, rad, mm/ms, rad/ms]
    end
    return y2, T, τfun, τifun
end

function cu.paramLumped(m::Wing2DOFModel, param::AbstractArray)
    cbar2, τ1, mwing, wΨ, τ2, Aw, dt = param
    kΨ, bΨ = hingeParams(wΨ)
    AR = Aw/cbar2
    cbar = sqrt(cbar2)
    Lw = Aw/cbar
    ycp = Lw*m.r1h # approx as in [Chen (2016)]
    return [one(eltype(param)), kΨ, bΨ, Aw^2*AR, mwing*ycp^2, mwing*cbar*ycp, mwing*cbar2], τ1, τ2, dt
end

function cu.paramAffine(m::Wing2DOFModel, opt::cu.OptOptions, traj::AbstractArray, param::AbstractArray, POPTS::cu.ParamOptOpts, scaleTraj=1.0; debugComponents::Bool=false)
    ny, nu, N, δt, liy, liu = cu.modelInfo(m, opt, traj)

    yo = k -> traj[liy[:,k]]*scaleTraj # traj is in output coords already
    uk = k -> @view traj[liu[:,k]]

    # Param stuff
    cbar2, τ1, mwing, wΨ, τ2, Aw, dtold = param
    AR = Aw/cbar2

    # THESE FUNCTIONS USE OUTPUT COORDS -------------
    # Nondimensionalize (wrt time) the velocities https://github.com/avikde/robobee3d/pull/119#issuecomment-577350049
    """Inertial"""
    # M = [Izz + mwing*ycp^2 + 1/2*cbar^2*mwing*γ^2*(1-cos(2*Ψ)) + m.ma/T^2   γ*cbar*mwing*ycp*cos(Ψ);  
        # γ*cbar*mwing*ycp*cos(Ψ)     Izz+cbar^2*γ^2*mwing]
    # Implicitly use Izz=0, Ixx=cbar^2 mw gamma^2
    function HMqTWithoutCoupling(ypos, yvel)
        φ, Ψ = ypos[1:2]
        dφ, Ψ̇ = yvel[3:4] * dtold
        return [0   0   0   0   dφ   0  1/2*γ^2*(1-cos(2*Ψ))*dφ     dφ*m.ma   dφ*m.ma*(-φ^2);
                0   0   0   0   0   0   2*Ψ̇*γ^2    0   0]
    end
    function HMqTCoupling(ypos, yvel)
        φ, Ψ = ypos[1:2]
        dφ, Ψ̇ = yvel[3:4] * dtold
        return [0   0   0   0   0   γ*Ψ̇*cos(Ψ)   0   0   0;
        0   0   0   0   0   γ*dφ*cos(Ψ)   0    0   0]
    end
    HMqT(ypos, yvel) = HMqTWithoutCoupling(ypos, yvel) + HMqTCoupling(ypos, yvel)

    """Coriolis"""
    function HC(y)
        # cor1 = [cbar^2*mwing*γ^2*sin(2*Ψ)*dφ*Ψ̇ - γ*cbar*mwing*ycp*sin(Ψ)*Ψ̇^2, 
        #     -cbar^2*mwing*γ^2*cos(Ψ)*sin(Ψ)*dφ^2]
        φ, Ψ = y[1:2]
        dφ, Ψ̇ = y[3:4] * dtold
        return m.bCoriolis ? [0   0   0   0   0   -γ*Ψ̇^2*sin(Ψ)   γ^2*sin(2*Ψ)*dφ*Ψ̇    0   0;
        0   0   0   0   0   0   -γ^2*cos(Ψ)*sin(Ψ)*dφ^2   0   0]  : zeros(2,9)
    end
    """Stiffness output"""
    function Hg(y)
        # [(m.ko*φ + m.ka/T*τifun(φ)), kΨ*Ψ]
        # τdamp = [-(m.bo + m.ba/T^2) * dφ, -bΨ * Ψ̇]
        φ, Ψ = y[1:2]
        return [m.kbo[1]*φ   0   0   0   0   0   0    0   0;
        0   Ψ   0    0   0   0   0   0   0]
    end
    """Stiffness actuator"""
    function Hgact(y)
        φ, Ψ = y[1:2]
        return [0   0   0   0   0   0   0    m.ka*φ   m.ka*(-φ^3/3);
        0   0   0   0   0   0   0   0   0]
    end
    """Damping"""
    function Hdamp(y)
        φ, Ψ = y[1:2]
        dφ, Ψ̇ = y[3:4] * dtold
        return [m.kbo[2]*dφ   0   0   0   0   0   0    m.ba*dφ   m.ba*dφ*(-φ^2);
        0   0   Ψ̇   0   0   0   0   0   0]
    end
    """Ext force (aero)"""
    function HF(y, tauaero)
        # Only keeping Fext_pdep=true version, and approximating only Faero propto this (ignore dependence of COP moving)
        tautil = -tauaero*dtold^2/(Aw^2*AR)

        # 1st order version
        return [0   0   0   tautil[1]   0   0   0   0   0;
        0   0   0   tautil[2]   0   0   0   0   0]
    end
    function HF2(y)
        _, Jaero, Faero = w2daero(m, y, param)
        return HF(y, Jaero'*Faero)
    end
    function Hvel(y)
        # Such that Hvel*pt[middle i.e./dt] = act. frame vel.
        φ, Ψ = y[1:2]
        dφ, Ψ̇ = y[3:4] * dtold
        return [0   0   0   0   0   0   0    dφ   dφ*(-φ^2)]
    end

    if debugComponents
        return yo, HMqTWithoutCoupling, HMqTCoupling, HC, Hg, Hgact, HF2, Hdamp, Hvel
    end
    
    # Functions to output
    "Takes in a Δy in output coords"
    function Hk(k, Δyprev, Δyk, Δykp1)
        y = yo(k) + Δyk
        ynext = yo(k+1) + Δykp1
        yprev = POPTS.centralDiff ? yo(max(k-1,1)) + Δyprev : y # Δyprev argument is ignored (and does not appear in jacobian)
        # Same Faero as before?
        _, Jaero, Faero = w2daero(m, y, param)
        # NOTE: it uses param *only for Faero*. Add same F as the traj produced
        # Divide up by dt-dependence https://github.com/avikde/robobee3d/pull/119#issuecomment-577350049
        H_dt2 = (HMqT(y, ynext) - HMqT(y, yprev))/(POPTS.centralDiff ? 2 : 1) + HC(y) + HF(y, Jaero'*Faero)
        H_dt1 = Hdamp(y)
        H_dt0 = Hg(y) + Hgact(y)
        # With nonlinear transmission need to break apart H
        φo = (yo(k) + Δyk)[1]
        return [cu.Hτ(H_dt2, φo)  cu.Hτ(H_dt1, φo)   cu.Hτ(H_dt0, φo)], cu.Hτ(Hvel(y), φo)
    end
    
    # For a traj, H(yk, ykp1, Fk) * pb = B uk for each k
    B = reshape([1.0, 0.0], 2, 1)

    return Hk, yo, uk, B, N
end

avgLift(m::Wing2DOFModel, opt::cu.OptOptions, traj::AbstractArray, param::AbstractArray) = mean(filter(!isnan, trajAero(m, opt, traj, param, :lift))) # last one is NaN for plotting
