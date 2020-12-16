
dt = 1/5000;

% plotTrial('../../../../Desktop/umpc_mm_data/trial36.mat', 1.5, dt)
% plotTrial('../../../../Desktop/umpc_mm_data/trial37.mat', 2, dt)
% plotTrial('../../../../Desktop/umpc_mm_data/trial32.mat', 1.9, dt)
% plotTrial('../../../../Desktop/umpc_mm_data/trial31.mat', 1.4, dt)
% plotTrial('../../../../Desktop/umpc_mm_data/trial21.mat', 1.5, dt)
plotTrial('../../../../Desktop/umpc_mm_data/trial20.mat', 1.4, dt)

function plotTrial(fname, tmax, sampling_time)
	load(fname);
	quickPlots(yout, sampling_time, tmax)
end
