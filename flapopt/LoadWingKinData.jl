using MAT, DSP, Dierckx
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


function analyzeData(fname)
	data, currTest = filterDataFull(fname, 2.2, 2.3, 600)
	# These are the channels
	drag, lift, opDisp, powerBias, powerSig, camTrig, sig, bias = data

	ld = plot(lift[:,1], lift[:,2], label="lift")
	plot!(ld, drag[:,1], -drag[:,2], label="drag")

	u = plot(sig[:,1], sig[:,2], label="sig")
	plot!(u, bias[:,1], bias[:,2], label="bias")

	pow = plot(powerSig[:,1], powerSig[:,2], label="psig")
	plot!(pow, powerBias[:,1], powerBias[:,2], label="pbias")

	p4 = plot(opDisp[:,1], opDisp[:,2], label="opDisp")
	# plot!(p4, camTrig[:,1], camTrig[:,2], label="camTrig")

	plot(ld, u, pow, p4)
end

analyzeData("../../../Desktop/vary_amplitude_no_lateral_wind_data/Test 22, 02-Sep-2016-11-39.mat")
