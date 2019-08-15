
include("OptBase.jl") #< including this helps vscode reference the functions in there

#=========================================================================
Param opt
=========================================================================#

function paramopt(m::Model, opt::OptOptions, traj::AbstractArray, params0::AbstractArray, εs; μs::Array{Float64}=[1e-1], Ninner::Int=1)
	ny, nu, N, δt, liy, liu = modelInfo(m, opt, traj)
	wkp = OptWorkspace(length(params0), (N+2)*ny)
	
	# FIXME: this is really more suited to the custom solver where Dg is already computed
	g_L, g_U = gbounds(m, opt, traj, εs...)
	Nc = length(g_L)
	Nx = length(traj)
	Dg = zeros(Nc, Nx)
	# Get Dg using the IPOPT functions
	nnz = Dgnnz(m, opt, traj)
	row = zeros(Int32, nnz)
	col = zeros(Int32, nnz)
	value = zeros(nnz)
	Dgsparse!(row, col, value, m, opt, traj, params0, :Structure, ny, nu, N, δt)
	Dgsparse!(row, col, value, m, opt, traj, params0, :Values, ny, nu, N, δt)
	for i = 1:nnz
		# DUMB
		Dg[row[i], col[i]] = value[i]
	end

	# Desired δx: lower u used
	δx = copy(traj)
	fill!(δx[1:(N+1)*ny], 0.0)
	δx[(N+1)*ny+1:(N+1)*ny+N*nu] .= -δx[(N+1)*ny+1:(N+1)*ny+N*nu] # negative of the currently applied force

	# Gradient wrt params TODO: more insight than autograd?
	function gofp(pp)
		gg = Array{Any,1}(undef, length(g_L))
		gvalues!(gg, m, opt, traj, pp, traj[1:ny])
		return gg
	end
	dg_dp = ForwardDiff.jacobian(gofp, params0)

	# Non-QP version first
	δp = -dg_dp \ Dg * δx
	println(δp)
end
