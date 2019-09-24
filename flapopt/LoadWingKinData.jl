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
	# vars["data"]

	# TODO:
	
	return vars["data"], vars["currTest"]
end


function analyzeData(fname)
	data, currTest = filterDataFull(fname, 2.2, 3.5, 600)

	t_lift = data[:,3]
	lift = data[:,4]

	t_drag = data[:,1]
	drag = data[:,2]

	t_opDisp = data[:,5]
	opDisp = data[:,6]

	t_powerBias = data[:,7]
	powerBias = data[:,8]

	t_powerSig = data[:,9]
	powerSig = data[:,10]
	t_camTrig = data[:,11]
	camTrig = data[:,12]

	t_sig = data[:,13]
	sig = data[:,14]

	t_bias = data[:,15]
	bias = data[:,16]

	pl1 = plot(t_lift, lift)
	plot(pl1)
end

analyzeData("../../../Desktop/vary_amplitude_no_lateral_wind_data/Test 22, 02-Sep-2016-11-39.mat")
