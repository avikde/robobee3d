using MAT
using Plots; gr()

"Has data for 8 channels, in pairs of (ti,datai) and the columns are stacked"
function filterDataFull(fname, tstart, tend, cutoff_freq)
	file = matopen(fname)
	vars = read(file)
	# display(vars["currTest"])
	# display(vars["currTest"]["Actuators"])
	sampleRate = vars["currTest"]["SampleRate"]
	freq = vars["currTest"]["Actuators"]["Frequency"]
	data = vars["data"]

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

		# return Nx2
		return [time_c sig_c]
	end

	dataOut = [filterChannel(i) for i = 1:8]
	
	return dataOut, vars["currTest"]
end


function analyzeData(fname)
	data, currTest = filterDataFull(fname, 2.2, 3.5, 600)
	# These are the channels
	drag, lift, opDisp, powerBias, powerSig, camTrig, sig, bias = data

	# t_lift = data[:,3]
	# lift = data[:,4]

	# t_drag = data[:,1]
	# drag = data[:,2]

	# t_opDisp = data[:,5]
	# opDisp = data[:,6]

	# t_powerBias = data[:,7]
	# powerBias = data[:,8]

	# t_powerSig = data[:,9]
	# powerSig = data[:,10]
	# t_camTrig = data[:,11]
	# camTrig = data[:,12]

	# t_sig = data[:,13]
	# sig = data[:,14]

	# t_bias = data[:,15]
	# bias = data[:,16]

	pl1 = plot(lift[:,1], lift[:,2])
	plot(pl1)
end

analyzeData("../../../Desktop/vary_amplitude_no_lateral_wind_data/Test 22, 02-Sep-2016-11-39.mat")
