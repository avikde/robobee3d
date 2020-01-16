using DelimitedFiles, Plots, LinearAlgebra, Statistics

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

alldata = vcat(transmissionStaticProcessCSV("poke1.csv"), transmissionStaticProcessCSV("poke2.csv"), transmissionStaticProcessCSV("poke3.csv"), transmissionStaticProcessCSV("poke4.csv"))
alldata = dataFilter(alldata)

mx, my, stdy, ycenter = discretize(alldata[:,1], alldata[:,2])
alldata[:,2] .-= ycenter

pl1 = scatter(alldata[:,1], alldata[:,2], xlabel="act", ylabel="tau")
plot!(pl1, mx, my, lw=2)

println(my)

# xs = range(1,stop=100, length=50)
# μs = log.(xs)
# σs = rand(length(xs))

pl2 = plot(mx, my, grid=true, ribbon=stdy, fillalpha=.5)

plot(pl1, pl2)
