
set(0, 'DefaultLineLineWidth', 1.5)

figure
clf
ax1 = gca();

figure('position', [0,0,600,250])
clf
ax2 = subplot(1,2,1);
ax3 = subplot(1,2,2);

plotTrial(ax1, ax2, ax3, '../../../../Desktop/umpc_mm_data/trial36.mat', 1.5, [0,0,1])
hold(ax1, 'all')
plotTrial(ax1, ax2, ax3, '../../../../Desktop/umpc_mm_data/trial37.mat', 1.6, [1,0,0]) % 2
plotTrial(ax1, ax2, ax3, '../../../../Desktop/umpc_mm_data/trial32.mat', 1.6, [0,0.5,0]) % 1.9
plotTrial(ax1, ax2, ax3, '../../../../Desktop/umpc_mm_data/trial31.mat', 1.4, [0.5,0.5,0])
plotTrial(ax1, ax2, ax3, '../../../../Desktop/umpc_mm_data/trial21.mat', 1.5, [0,0.5,0.5])
plotTrial(ax1, ax2, ax3, '../../../../Desktop/umpc_mm_data/trial20.mat', 1.4, [0.5,0,0.5])

hold(ax1, 'off')
grid(ax1, 'on')
view(ax1, [-35,8])
pbaspect(ax1,[1,1,1])

xlabel(ax1, 'x [m]')
ylabel(ax1, 'y [m]')
zlabel(ax1, 'z [m]')

plot(ax2, [0,2], [0,0], 'k--')
ylim(ax2, [-pi/4, pi/4])
plot(ax3, [0,2], [0,0], 'k--')
ylim(ax3, [-pi/4, pi/4])
ylabel(ax2, 'Roll [rad]')
ylabel(ax3, 'Pitch [rad]')
xlabel(ax2, 'Time [s]')
xlabel(ax3, 'Time [s]')

function plotTrial(ax1, ax2, ax3, fname, tmax, col)
	%dt = 1/5000;
	load(fname);
	tt = yout(:,1);
	qvicon = yout(:,7:12);
	plog = yout(:,16:18);
	dqlog = yout(:,19:24);
	uquadlog = yout(:,26:28);
	accdeslog = yout(:,29:34);
	%quickPlots(yout, sampling_time, tmax)
	[s, ds] = uprightCalcs(qvicon(:,4:6), dqlog(:,4:6));
	
	traj3plot(ax1, ax2, ax3, tt, qvicon, s, tmax, col)
end
