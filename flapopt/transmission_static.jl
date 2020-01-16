using DelimitedFiles, Plots, LinearAlgebra, Statistics, Dierckx

fname = "poke1.csv"

function transmissionStaticProcessCSV(fname)
	dat = readdlm(string("data/transmission/",fname), ',', Float64, skipstart=2)
	# cols are t, Ax,Ay, Bx,By, Cx,Cy
	# Ax is the x pos of the actuator [mm]
	# B, C are points on the wing (arctan should be wing angle)

	dat[:,2] .-= dat[1,2]
	Np = size(dat, 1)
	wingAngs = zeros(Np)
	for i=1:Np
		wingVec = dat[i,6:7] - dat[i,4:5]
		wingAngs[i] = atan(wingVec[2]/wingVec[1])
	end

	return hcat(dat[:,2], wingAngs)
end

function dataFilter(alldata)
	# points on the left extreme edge (also makes the xdata symmetric)
	Np = size(alldata, 1)
	return hcat([alldata[i,:] for i=1:Np if abs(alldata[i,1]) < 0.2]...)'
end

function discretize(xdata, ydata, nbins=30)
	binedges = range(minimum(xdata), stop=maximum(xdata), length=nbins+1)
	Np = length(xdata)
	# arrays to store the mean, stddev
	xs = [Float64[] for i=1:nbins]

	for i in 1:nbins
		for j = 1:Np
			if xdata[j] >= binedges[i] && xdata[j] < binedges[i+1]
				# accumulate the data and the number of data pts in each bin
				xs[i] = [xs[i]; ydata[j]]
			end
		end
	end

	# Statistics
	means = [mean(xs[i]) for i=1:nbins]
	stds = [std(xs[i]) for i=1:nbins]
	bincenters = 0.5 * (binedges[1:end-1] + binedges[2:end])

	# y data mean should be zero (video tilted) - shift by y of center binedges
	ycenter = 0.5*(means[nbins÷2] + means[nbins÷2+1])
	means .-= ycenter

	return bincenters, means, stds, ycenter
end

function predictedTransmission(tau1, tau2, ycp)
	# tau1,tau2 for the 2D model; need to divide by ycp
	tau(sigmaa) = (tau1*sigmaa + tau2*sigmaa^3/3)/ycp
	Dtau(sigmaa) = (tau1 + tau2*sigmaa^2)/ycp

	return tau, Dtau
end

function measDtau(x, fx)
	# DSP.jl
	# myfilter = digitalfilter(Lowpass(1), Butterworth(10))
	# fy = filtfilt(myfilter, y)
	# A_x = 1.:2.:40.
	# A = [log(x) for x in A_x]

	# # Interpolations.jl
	# itp = interpolate(fx, BSpline(Cubic(Line(OnGrid()))))
	# sitp = scale(itp, x)
	# g = vcat([Interpolations.gradient(sitp, xi) for xi in x]...)
	# return sitp.(x), g

	# using Dierckx.jl
	spl = Spline1D(x, fx)
	Df = [derivative(spl, xi) for xi in x]
	
	return spl.(x), Df
end

# Predicted transmission functions
tau1 = 15.98
tau2 = 31.959
ycp = 9.1
tau, Dtau = predictedTransmission(tau1, tau2, ycp)
taul, Dtaul = predictedTransmission(tau1, 0, ycp)

# Load data
alldata = vcat(transmissionStaticProcessCSV("poke1.csv"), transmissionStaticProcessCSV("poke2.csv"), transmissionStaticProcessCSV("poke3.csv"), transmissionStaticProcessCSV("poke4.csv"))
alldata = dataFilter(alldata)

mx, my, stdy, ycenter = discretize(alldata[:,1], alldata[:,2])
alldata[:,2] .-= ycenter
itpy, itpDy = measDtau(mx, my)

# Scatter plot of data
pl1 = scatter(alldata[:,1], alldata[:,2], xlabel="act", ylabel="tau")
plot!(pl1, mx, my, lw=2)

# println(my)

# Overlaid with predicted
xx = -0.3:0.01:0.3
pl2 = plot(xx, -tau.(xx), grid=true, lw=2, title="Transmission function", xlabel="act disp [mm]", ylabel="output [rad]", label="Nonlin")
plot!(pl2, xx, -taul.(xx), lw=2, label="Lin")
plot!(pl2, mx, my, lw=2, ribbon=stdy, fillalpha=.3, label="Meas")
plot!(pl2, mx, itpy, lw=2, label="itp")

pl3 = plot(xx, Dtau.(xx), lw=2, title="Instantaneous tr. ratio", xlabel="act disp [mm]", ylabel="dwing/dact [rad/mm]", label="Nonlin")
plot!(pl3, xx, Dtaul.(xx), lw=2, label="Lin")
plot!(pl3, mx, -itpDy, lw=2, label="itp")

plot(pl2, pl3)
