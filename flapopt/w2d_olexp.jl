using Plots, DelimitedFiles

mod4bh1_V = [#= 150, =# 160, 170, 180, 190, 200]

mod4bh1 = [
	# [150	47.08;
	# 160	61.52;
	# 170	75.67;
	# 175	75.99;
	# 180	73.98],
	[150	52.05;
	160	72.76;
	170	81.88;
	175	80.78;
	180	76.77],
	[150	57.47;
	160	82.13;
	165	84.7;
	170	86.7;
	180	78.83],
	[140	43.96;
	150	58.41;
	160	72.68;
	# 170	64.02;
	175	86.21;
	180	82.02;
	185	79.17;
	190	73.23],
	[140	46.51;
	150	72.35;
	155	87.65;
	158	92.36;
	160	93.17;
	165	65.21;
	170	70.31],
	[140	47.67;
	150	73.97;
	155	90.75;
	157	95.66;
	160	82.74;
	170	74.94]
	]

function readOLExpCSV(fname)
	dat = readdlm(fname, ',', Float64, skipstart=1)
	Vs = unique(dat[:,1])
	return Vs, [dat[dat[:,1] .== V,2:3] for V in Vs]
end

function olExpPlot(V, stroke; showV=[])
	Nv = length(V)
	p = plot(xlabel="Freq [Hz]", ylabel="Norm. stroke ampl [deg/V]", ylims=(0.2,0.7))
	for i=1:Nv
		if length(showV) == 0 || Int(V[i]) in showV
			plot!(p, stroke[i][:,1], stroke[i][:,2]/V[i], label=Int(V[i]), markershape=:auto, lw=2)
		end
	end
	return p
end

mod1a1_V, mod1a1 = readOLExpCSV("data/normstroke/Param opt manuf 2 - mod1 a1 redo.csv")

plot(
	olExpPlot(mod1a1_V, mod1a1; showV=[160,170,180,190,200]), 
	olExpPlot(mod4bh1_V, mod4bh1))
gui()
	