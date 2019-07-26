
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
struct Wing2DOFModel <: controlutils.Model end

# Fixed params -----------------
const R = 20


function cu.dims(m::Wing2DOFModel)::Tuple{Int, Int}
    return 4, 1
end

function cu.limits(m::Wing2DOFModel)::Tuple{Vector, Vector, Vector, Vector}
    # This is based on observing the OL trajectory. See note on units above.
    umax = @SVector [75] # [mN]
    umin = -umax
    xmax = @SVector [300e-3, 1.5, Inf, Inf] # [mm, rad, mm/ms, rad/ms]
    xmin = -xmax
    return umin, umax, xmin, xmax
end

function cu.limitsTimestep(m::Wing2DOFModel)::Tuple{Float64, Float64}
	return 0.01, 10.0
end


"Returns paero [mm], Jaero, Faero [mN]"
function w2daero(y::Vector, u::Vector, _params::Vector)
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
    paero = wing1 + RΨ * @SVector [0, -cbar]
    # TODO: add moving CoP from [Chen et al (2017)]
    Jaero = @SMatrix [1 cbar * cΨ; 0 cbar * sΨ]
    
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
    α = π/2 - Ψ # AoA
    Caero = @SVector [((CDmax + CD0)/2 - (CDmax - CD0)/2 * cos(2α)), CLmax * sin(2α)]
    Faero = 1/2 * ρ * cbar * R * σ̇^2 * Caero * sign(-σ̇) # [mN]

    return paero, Jaero, Faero
end

"Continuous dynamics second order model"
function cu.dydt(model::Wing2DOFModel, y::Vector, u::Vector, _params::Vector)::Vector
    # unpack
    cbar, T = _params
    σ, Ψ, σ̇, Ψ̇ = (@SVector [T, 1, T, 1]) .* y # [mm, rad, mm/ms, rad/ms]
    # NOTE: for optimizing transmission ratio
    # Thinking of y = (sigma_actuator, psi, dsigma_actuator, dpsi)
    # u = (tau_actuator)
    # sigma = sigma_actuator * T; tau = tau_actuator / T
    cΨ = cos(Ψ)
    sΨ = sin(Ψ)

    # params
    mspar = 0 # [mg]
    mwing = 0.51 # [mg]
    Iwing = mwing * cbar^2 # cbar is in mm
    kσ = 0
    kΨ = 5 # [mN-mm/rad]
    bΨ = 1 # [mN-mm/(rad/ms)]

    # inertial terms
    M = @SMatrix [mspar+mwing   cbar*mwing*cΨ; cbar*mwing*cΨ   Iwing+cbar^2*mwing]
    corgrav = @SVector [kσ*σ - cbar*mwing*sΨ*Ψ̇^2, kΨ*Ψ]
    # non-lagrangian terms
    τdamp = @SVector [0, -bΨ * Ψ̇]
    _, Jaero, Faero = w2daero(y, u, _params)
    τaero = Jaero' * Faero # units of [mN, mN-mm]
    # input
    τinp = @SVector [u[1]/T, 0]

    ddq = inv(M) * (-corgrav + τdamp + τaero + τinp)
    # return ddq
    return [y[3], y[4], ddq[1], ddq[2]]
end

# Create an initial traj --------------

"""
freq [kHz]; posGains [mN/mm, mN/(mm-ms)]; [mm, 1]
Example: trajt, traj0 = Wing2DOF.createInitialTraj(0.15, [1e3, 1e2], params0)
"""
function createInitialTraj(m::Wing2DOFModel, N::Int, freq::Real, posGains::Vector, params0::Vector)
    # Create a traj
    σmax = cu.limits(m)[end][1]
    function strokePosController(y, t)
        σdes = 0.9 * σmax * sin(freq * 2 * π * t)
        return posGains[1] * (σdes - y[1]) - posGains[2] * y[3]
    end
    strokePosControlVF(y, p, t) = cu.dydt(m, y, [strokePosController(y, t)], params0)
    # OL traj1
    teval = collect(0:1e-1:100) # [ms]
    prob = ODEProblem(strokePosControlVF, zeros(4), (teval[1], teval[end]))
    sol = solve(prob, saveat=teval)

    # σt = plot(sol, vars=3, ylabel="act vel [m/s]")
    # Ψt = plot(sol, vars=2, ylabel="hinge ang [r]")
    # plot(σt, Ψt, layout=(2,1))
    # gui()

    starti = 172
    olRange = starti:3:(starti + 3*N)
    trajt = sol.t[olRange]
    δt = trajt[2] - trajt[1]
    olTrajaa = sol.u[olRange] # 23-element Array{Array{Float64,1},1} (array of arrays)
    olTraju = [strokePosController(olTrajaa[i], trajt[i]) for i in 1:N] # get u1,...,uN
    traj0 = [vcat(olTrajaa...); olTraju; δt] # dirtran form {x1,..,x(N+1),u1,...,u(N),δt}

    return trajt .- trajt[1], traj0
end

function cu.plotTrajs(m::Wing2DOFModel, t::Vector, params::Vector, args...; vart=true)
	ny, nu = cu.dims(m)
	traj = args[1]

	N = cu.Nknot(m, traj; vart=vart)
    liy, liu = cu.linind(m, N)
	Ny = (N+1)*ny
	# stroke "angle" = T*y[1] / R
	cbar, T = params
	σt = plot(t, traj[@view liy[1,:]] * T / (R/2), marker=:auto, ylabel="stroke ang [r]", legend=false, title="timestep=$(round(traj[end]; sigdigits=4))ms; c=$(round(cbar; sigdigits=4))mm, T=$(round(T; sigdigits=4))")
	Ψt = plot(t, traj[@view liy[2,:]], marker=:auto, legend=false, ylabel="hinge ang [r]")
    ut = plot(t, [traj[@view liu[1,:]];NaN], marker=:auto, legend=false, ylabel="stroke force [mN]")
    # Plot of aero forces at each instant
    Faerok = k -> w2daero(traj0[@view liy[:,k]], traj0[@view liu[:,k]], params0)[end]
    Faeros = hcat([Faerok(k) for k=1:N]...)
	aerot = plot(t, [Faeros';[NaN NaN]], marker=:auto, legend=false, ylabel="aero force [mN]")
	plot(σt, Ψt, ut, aerot)
	gui()
end

# "Cost function components" ------------------

"Objective to minimize"
function cu.Jobj(m::Wing2DOFModel, traj::Vector, params::Vector; vart::Bool=true, fixedδt::Float64=1e-3)::Number
	N = cu.Nknot(m, traj; vart=vart)
    liy, liu = cu.linind(m, N)
	δt = vart ? traj[end] : fixedδt
    Favg = @SVector zeros(2)
    for k = 1:N
        vy = @view liy[:,k]
        vu = @view liu[:,k]
        paero, _, Faero = w2daero(traj[vy], traj[vu], params)
        Favg += Faero
    end
    # Divide by the total time
    Favg /= (N*δt) # [mN]
    # max avg lift
	return -Favg[2]
end


