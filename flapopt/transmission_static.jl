using DelimitedFiles, Plots

fname = "data/transmission/poke1.csv"
dat = readdlm(fname, ',', Float64, skipstart=2)
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

pl1 = scatter(dat[:,2], wingAngs, xlabel="act", ylabel="wing")
plot(pl1)
