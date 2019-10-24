using MAT, DSP, Dierckx, DelimitedFiles, LinearAlgebra, Statistics
using Plots; gr()

"Has data for 8 channels, in pairs of (ti,datai) and the columns are stacked"
function filterDataFull(fname, tstart, tend, cutoff_freq)
	file = matopen(fname)
	vars = read(file)
	# display(vars["currTest"])
	# display(vars["currTest"]["Actuators"])
	sampleRate = vars["currTest"]["SampleRate"]
	freq = vars["currTest"]["Actuators"]["Frequency"][1]
	data = vars["data"]

	# Design the filter
	myfilter = digitalfilter(Lowpass(cutoff_freq; fs=sampleRate), Butterworth(10))

	function filterChannel(i)
		# shorten length of 8 channels to eliminate ramp up/ ramp down 
		time_ci = data[:,2*i-1] #time_ci is time of channel i
		sig_c = data[:,2*i] #sig_ci is signal of channel i
		ind_tstart = findfirst(x -> x >= tstart, time_ci)
		ind_tend = findfirst(x -> x > tend, time_ci)
		
		time_c = time_ci[ind_tstart:ind_tend]
		sig_c = sig_c[ind_tstart:ind_tend]
		
		# 	%specifically for lift, subtract off mean
		# 	if i==1 || i==2
		# 		runTime =5;
		# 		len =currTest.SampleRate*runTime;
		# 		%mean([data(1:round(len/20),2); data(round(19*len/20):len,2)])
		# 		sig_c = -0.912*(sig_c - mean([data(1:round(len/20),2*i); data(round(19*len/20):len,2*i)]))*1000; %convert to unit of mg 
			
		# %        plot(time_c,sig_c,'r-')
		# %        plot(data(:,2*i-1),data(:,2*i),'g-')
		# %        p=5

		sig_c_mod = filtfilt(myfilter, sig_c)
		
		dt = 1/(freq*30) # 30 points per cycle

		t_start = time_c[1]+5/freq # throw away first five cycles
		cycles = floor((time_c[end]-time_c[1])*freq) - 5 # throw away last five cycles
		t_end = t_start + cycles * (1/freq)

		# return Nx2
		spl = Spline1D(time_c, sig_c_mod; k=2)
		tt = t_start:dt:t_end
		return [tt   spl(tt)]
	end

	dataOut = [filterChannel(i) for i = 1:8]
	
	return dataOut, vars["currTest"]
end


function loadDAQData(fname)
	data, currTest = filterDataFull(fname, 2.2, 2.3, 600)
	# These are the channels
	drag, lift, opDisp, powerBias, powerSig, camTrig, sig, bias = data

	# ld = plot(lift[:,1], lift[:,2], label="lift")
	# plot!(ld, drag[:,1], -drag[:,2], label="drag")

	# u = plot(sig[:,1], sig[:,2], label="sig")
	# plot!(u, bias[:,1], bias[:,2], label="bias")

	# pow = plot(powerSig[:,1], powerSig[:,2], label="psig")
	# plot!(pow, powerBias[:,1], powerBias[:,2], label="pbias")

	# p4 = plot(opDisp[:,1], opDisp[:,2], label="opDisp")
	# plot!(p4, camTrig[:,1], camTrig[:,2], label="camTrig")

	# plot(ld, u, pow, p4)

	return Dict("drag" => drag, "lift" => lift, "sig" => sig), currTest
end

"""Align the rows at the bottom so the all have the same #rows first.
Assumes the first 2 rows are text and headers.
Points A, B should be 2 distinct points on the spar. Points C, D can be different points not on the spar.
dC, dD = (perp) distances of those points from the wing spar (this is the only metric information needed).
vidX = vid sample rate / 30 (or whatever was chosen when exporting avi). e.g. 7500fps -> 30 => vidX = 250
trialFreq = freq of flapping (used *only* to trim to 1 flap)
"""
function loadVideoData(fname; dC=5.345, dD=4.047, vidX=250, trialFreq=165, makegif=false)
	dat = readdlm(fname, ',', Float64, skipstart=2)
	# Should be an Nx12 array, for mass A (t, x, y), ... mass D
	# Use the first col as the time vector
	tq = dat[:,1]

	"i = massid between 1 and 4"
	function xy_at_t(i; kwargs...)
		splx = Spline1D(dat[:,3*i-2], dat[:,3*i-1]; kwargs...)
		sply = Spline1D(dat[:,3*i-2], dat[:,3*i]; kwargs...)
		return [splx(tq) sply(tq)]
	end

	"Find p0, the center of rotation assuming pA, pB in the same plane using linear LS"
	function find_p0(pA, pB)
		A = zeros(Np*2,2)
		b = zeros(Np*2)
		for i in 1:Np
			# these are row vectors
			_pA = pA[i,:]
			_pB = pB[i,:]
			pABperp = [0 -1; 1 0] * (_pB - _pA)
			# Append to the A,b mats for LS
			A[2*i-1:2*i,:] = [pABperp'; pABperp']
			b[2*i-1:2*i] = [dot(pABperp, _pA); dot(pABperp, _pB)]
		end

		# solve A*p0 = b
		return A\b
	end

	function find_Φ(pAi, pBi, p0)
		A = [(pAi - p0)'; (pBi - p0)']
		# Want a nullspace vector of A https://cseweb.ucsd.edu/classes/wi15/cse252B-a/nullspace.pdf
		F = svd(A)
		# Last singular value in F.S should be small. Corresponding right nullsp vector is the last col of V
		cΦ, sΦ = F.Vt[end,:]
		# This option likely depends on whether a coordinate frame was introduced in Tracker.
		# return atan(sΦ/cΦ)-π/2
		return atan(-cΦ/sΦ) + π # no coord frame in tracker (pixel coords)
	end

	function find_Ψ(pBi, pCi, pDi, p0, Φi)
		# First some helper functions
		Rot = Φ -> [cos(Φ) -sin(Φ); sin(Φ) cos(Φ)]

		function find_proj_dist(pOnSpar, pOnWing, p0, Φ)
			# project pC on to the spar line to get pC0
			pProjOnSpar = p0 + dot(pOnWing - p0, pOnSpar - p0) * (pOnSpar - p0) / dot(pOnSpar - p0, pOnSpar - p0)
			rotToSparFrame = Rot(-Φ) * (pOnWing - pProjOnSpar)
			return rotToSparFrame[2] # return the component perp to the spar line
		end

		c = find_proj_dist(pBi, pCi, p0, Φi)
		d = find_proj_dist(pBi, pDi, p0, Φi)

		# sine of the hinge angle appears in these projections
		sΨ = [dC; dD] \ [c;d]
		return asin(sΨ)
	end

	# get time in ms
	tms = 1000*tq/vidX
	pA, pB, pC, pD = [xy_at_t(i; k=1) for i=1:4] # Nx2 x 4
	Np = length(tq)
	p0 = find_p0(pA, pB)
	# println("p0 = ", p0)
	Φ = [find_Φ(pA[i,:], pB[i,:], p0) for i=1:Np]
	Ψ = [find_Ψ(pB[i,:], pC[i,:], pD[i,:], p0, Φ[i]) for i=1:Np]
	
	# Trim to tms=1000/trialFreq
	ind_1cyc = findfirst(x -> x >= 1000/trialFreq, tms)
	tms = tms[1:ind_1cyc]
	Φ = Φ[1:ind_1cyc]
	Ψ = Ψ[1:ind_1cyc]
	Np = length(tms)
	
	if makegif
		function drawFrame(k)
			span = norm(pB[1,:] - p0) * 1.5
			w = plot([p0[1]], [p0[2]], marker=:auto, color=:black, label="p0", xlims=(p0[1]-1.5*span, p0[1]+1.5*span), ylims=(p0[2]-1.5*span, p0[2]+1.5*span), aspect_ratio=1, legend=false) # unfortunately the legend screws up the aspect ratio if it is outer
			plot!(w, [pA[k,1]], [pA[k,2]], marker=:auto, color=:red, label="pA")
			plot!(w, [pB[k,1]], [pB[k,2]], marker=:auto, color=:cyan, label="pB")
			plot!(w, [pC[k,1]], [pC[k,2]], marker=:auto, color=:magenta, label="pC")
			plot!(w, [pD[k,1]], [pD[k,2]], marker=:auto, color=:purple, label="pD")
			# Stroke line
			plot!(w, [p0[1], p0[1] + span*cos(Φ[k])], [p0[2], p0[2] + span*sin(Φ[k])], color=:black, label="spar")

			w2 = plot(tms, Φ.-mean(Φ), linewidth=2, label="stroke")
			plot!(w2, tms, Ψ, linewidth=2, label="hinge")
			vline!(w2, [tms[k]])
			return plot(w, w2, layout=(2,1))
		end
		@gif for k = 1:Np
			drawFrame(k)
		end
	end

	return tms, Φ.-mean(Φ), Ψ
end

"""
vidX -- see loadAlignedData()
sigYout -- if true, assume yout from the target PC (assume 10KHz sample rate); if false, assume old Kevin data
ForcePerVolt -- actuator force mN/V

For sigYout = true:
- startOffs -- the first video frame in the CSV relative to the trigger (in PCC what you set as `[`)
- sigi -- which signal S1 or S2 (pass in 1 or 2)

For sigYout = false:
- startOffs -- the first timestamp used from the MAT file
"""
function loadAlignedData(fnameMat, fnameCSV, startOffs; strokeMult=1.0, ForcePerVolt=0.75, sigYout=true, vidSF=7500, sigi=1, sigsign=1.0)
	if sigYout
		# New yout created on the target computer. Assuming sample rate is 10KHz.
		# Column index
		# 1		Time
		# 2		Bias (ramp up ~ flat ~ ramp down)
		# 3 	Voltage_left
		# 4		Voltage_right
		# 5		Input voltage frequency
		# 6		dt (?)
		file = matopen(fnameMat)
		vars = read(file)
		yout = vars["yout"]
		freq = yout[1,5]

		# assuming t=0 coincides for both due to the trigger
		tstartMat = startOffs/vidSF # [s]
		sigt = yout[:,1]
		sig = yout[:,sigi+2] - 0.5 * yout[:,2] # subtract half bias

		# # Test plot
		# sigs = plot(yout[:,1], yout[:,2], label="b")
		# plot!(sigs, yout[:,1], yout[:,3], label="s1")
		# plot!(sigs, yout[:,1], yout[:,4], label="s2")
		# # others = plot(yout[:,1], yout[:,5], label="freq")
		# # plot!(others, yout[:,1], yout[:,6], label="col6")
		# plot(sigs)
		# gui()
	else
		# Old Kevin data
		daq, currTest = loadDAQData(fnameMat)
		freq = currTest["Actuators"]["Frequency"][1]
		Vpp = currTest["Actuators"]["Amplitude"][1]
		tstartMat = startOffs
		sigt = daq["sig"][:,1]
		sig = daq["sig"][:,2] .- Vpp * DAQ_To_Volts/2
	end
	tms, Φ, Ψ = loadVideoData(fnameCSV; trialFreq=freq, vidX=vidSF/30)

	# Align to get input --------------------------------------------
	
	# Find one cycle of data from the mat
	ind_tstart = findfirst(x -> x >= tstartMat, sigt)
	ind_tend = findfirst(x -> x >= tstartMat + 1/freq, sigt)

	function _alignDAQToVideo()
		# align voltage to the same times
		t = sigt[ind_tstart:ind_tend]
		t .= (t .- t[1]) * 1000 # to ms
		v = sig[ind_tstart:ind_tend]

		spl = Spline1D(t, v; k=2)
		return spl(tms)
	end
	# Vpp = max.(sig2)
	DAQ_To_Volts = 100 # Manually checked for 180V, max was 1.795
	uact = sigsign * _alignDAQToVideo() * DAQ_To_Volts * ForcePerVolt

	# Filter and take derivative of the video data -----------------------
	# sample rate
	fs = 1/mean(diff(tms))

	function numDeriv(v)
		dv = diff(v) ./ diff(tms)
		return [dv; dv[1]]
	end

	"cutoff_freq is in KHz"
	function lpfilt(v, ord, cutoff_freq)
		# Design the filter
		myfilter = digitalfilter(Lowpass(cutoff_freq; fs=fs), Butterworth(ord))
		return filtfilt(myfilter, v)
	end
	
	Ψ .= lpfilt(Ψ, 2, 1.5)
	dΦ = lpfilt(numDeriv(Φ), 2, 1.5)
	dΨ = lpfilt(numDeriv(Ψ), 2, 1.0)
	# ---------------------------------------------------------------------

	# Now convert to dirtran form for compatibilty with prior code
	Ndp1 = length(Φ)
	Y = [strokeMult*Φ';Ψ';strokeMult*dΦ';dΨ']
	X = [reshape(Y, 4*Ndp1);uact[1:end-1]]

	# aa = plot(tms, Φ)
	# plot!(aa, tms, Ψ)
	# bb = plot(tms, dΦ)
	# plot!(bb, tms, dΨ)
	# plot(aa, bb, layout=(2,1))
	# gui()
	return Ndp1-1, tms, X#, alignDAQToVideo(daq["lift"]), alignDAQToVideo(daq["drag"])
end

# loadAlignedData("../../../Desktop/vary_amplitude_no_lateral_wind_data/Test 22, 02-Sep-2016-11-39.mat", "data/lateral_windFri Sep 02 2016 18 45 18.344 193 utc.csv", 2.24)
