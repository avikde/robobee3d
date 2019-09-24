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

	return sig, currTest
end

"Align the rows at the bottom so the all have the same #rows first.
Assumes the first 2 rows are text and headers.
dC, dD are the distances of those points from the wing spar (this is the only metric information needed)."
function loadVideoData(fname; dC=1.0, dD=1.0, vidX=200, trialFreq=130)
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
		return atan(sΦ/cΦ)-π/2
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
	
	function drawFrame(k)
		span = 12.8
		w = plot([p0[1]], [p0[2]], marker=:auto, color=:black, label="p0", xlims=(0,30), ylims=(-5,10), aspect_ratio=1)
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
	# @gif for k = 1:Np
	# 	drawFrame(k)
	# end

	return tms, Φ.-mean(Φ), Ψ
end

"tstartMat is the first timestamp used from the MAT file, and 1 cycle is used"
function loadAlignedData(fnameMat, fnameCSV, tstartMat)
	sig, currTest = loadDAQData(fnameMat)
	freq = currTest["Actuators"]["Frequency"][1]
	Vpp = currTest["Actuators"]["Amplitude"][1]
	tms, Φ, Ψ = loadVideoData(fnameCSV; trialFreq=freq)
	
	# Find one cycle of data from the mat
	ind_tstart = findfirst(x -> x >= tstartMat, sig[:,1])
	ind_tend = findfirst(x -> x >= tstartMat + 1/freq, sig[:,1])

	function alignDAQToVideo(daqVec)
		# align voltage to the same times
		t = daqVec[ind_tstart:ind_tend,1]
		t .= (t .- t[1]) * 1000 # to ms
		v = daqVec[ind_tstart:ind_tend,2]

		spl = Spline1D(t, v; k=2)
		return spl(tms)
	end
	# Vpp = max.(sig2)
	DAQ_To_Volts = 100 # Manually checked for 180V, max was 1.795
	Volts_To_Force = 0.75 # [mN/V]
	uact = (alignDAQToVideo(sig) * DAQ_To_Volts .- Vpp/2) * Volts_To_Force
	return tms, Φ, Ψ, uact
end

# loadAlignedData("../../../Desktop/vary_amplitude_no_lateral_wind_data/Test 22, 02-Sep-2016-11-39.mat", "data/lateral_windFri Sep 02 2016 18 45 18.344 193 utc.csv", 2.24)
