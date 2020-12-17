
dt = 1/5000;

plotTrial('../../../../Desktop/umpc_mm_data/trial36.mat', 1.5, dt)
hold all
% plotTrial('../../../../Desktop/umpc_mm_data/trial37.mat', 2, dt)
% plotTrial('../../../../Desktop/umpc_mm_data/trial32.mat', 1.9, dt)
% plotTrial('../../../../Desktop/umpc_mm_data/trial31.mat', 1.4, dt)
% plotTrial('../../../../Desktop/umpc_mm_data/trial21.mat', 1.5, dt)
% plotTrial('../../../../Desktop/umpc_mm_data/trial20.mat', 1.4, dt)


hold off
grid on

function plotTrial(fname, tmax, sampling_time)
	load(fname);
	tt = yout(:,1);
	qvicon = yout(:,7:12);
	plog = yout(:,16:18);
	dqlog = yout(:,19:24);
	uquadlog = yout(:,26:28);
	accdeslog = yout(:,29:34);
	%quickPlots(yout, sampling_time, tmax)
	
	traj3plot(tt, qvicon(:,1:3), tmax)
end
