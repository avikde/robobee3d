using DelimitedFiles, Plots, LinearAlgebra

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

alldata = vcat(transmissionStaticProcessCSV("poke1.csv"), transmissionStaticProcessCSV("poke2.csv"), transmissionStaticProcessCSV("poke3.csv"), transmissionStaticProcessCSV("poke4.csv"))

pl1 = scatter(alldata[:,1], alldata[:,2], xlabel="act", ylabel="tau")

xs = range(1,stop=100, length=50)
μs = log.(xs)
σs = rand(length(xs))

pl2 = plot(xs,μs,grid=false,ribbon=σs,fillalpha=.5)

plot(pl1, pl2)
