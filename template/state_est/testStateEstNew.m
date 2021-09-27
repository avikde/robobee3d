% load('../../../../Desktop/20200926_closedloop/p10_d0_2.mat')
clear stateEst
clear
load('../../../../Desktop/uprightmpc2/open3.mat')
set(0, 'DefaultLineLineWidth', 2)
Np = size(yout, 1);
SAMPLE_RATE = 10000;
DEADTIME = 0.01;

% Run actual stateEst filter
dqfilt = zeros(Np, 6);

% 
s = zeros(Np,3);
sf = zeros(Np,3);
ds = zeros(Np,3);
dsf = zeros(Np,3);
dsnum = zeros(Np,3);
omgnum = zeros(Np,3); % from Rb WORKING
ww = zeros(Np, 3); % differentiate ZYZ angles NOT WORKING
e3h = [0 -1 0; 1 0 0; 0 0 0];
Rbprev = eye(3);
% Log parsing
qvicon = yout(:,7:12); %13--15 drive signals log
% State est outputs not really used
plog = yout(:,16:18); % state est outputs log
dqlog = yout(:,19:24);
% 25=drv_bias
% 26:28=uquad
% 29:34=accdes

for i=1:Np
% Run actual stateEst filter
	[p, Rbf, dq, nmeas] = stateEst(yout(i,1), yout(i,7:12)', DEADTIME);
	dqfilt(i,:) = dq;
	
	% Get Rb from raw vicon
	Rbvic = eul2rotm(qvicon(i,4:6), 'ZYZ');
	s(i,:) = Rbvic * [0 0 1]';
	sf(i,:) = Rbf * [0 0 1]';
	if i > 1
		dt = yout(i,1) - yout(i-1,1);
		sdiff = (s(i,:) - s(i-1,:)) / dt;
		dsnum(i,:) = dsnum(i-1,:) + 0.005 * (sdiff - dsnum(i-1,:));
		%
		Rbdot = (Rbvic - Rbprev) / dt;
		omghat = Rbvic' * Rbdot;
		omg1 = [omghat(3,2),omghat(1,3),omghat(2,1)];
		omgnum(i,:) = omgnum(i-1,:) + 0.005 * (omg1 - omgnum(i-1,:));
		%This is using state est outputs (not plotted)
		angdiff = (dqlog(i,4:6) - dqlog(i-1,4:6)) / dt;
		ww(i,:) = ww(i-1,:) + 0.005 * (angdiff - ww(i-1,:));
	end
	ds(i,:) = -Rbvic * e3h * omgnum(i,:)';
	dsf(i,:) = -Rbf * e3h * dqlog(i,4:6)';
	Rbprev = Rbvic;
end
fprintf(1, 'approx vicon rate discarding nan = %f\n', nmeas/Np * SAMPLE_RATE)

clf
subplot(221)
hold all
% plot(yout(:,1), dqfilt(:,4:6))
% hold all
plot(yout(:,1), qvicon(:,4:6))
%plot(yout(:,1), dqfilt(:,4:6), '--')
hold off
% legend('px','py','pz','fpx','fpy','fpz')

subplot(222)
hold all
plot(yout(:,1), omgnum)
%plot(yout(:,1), yout(:,25:27))
plot(yout(:,1), dqfilt(:,4:6), '--')
hold off
legend('wx','wy','wz','fwx','fwy','fwz')

subplot(223)
hold all
plot(yout(:,1), s)
plot(yout(:,1), sf, '--')% 
hold off
legend('sx','sy','sz','fsx','fsy','fsz')

subplot(224)
hold all
plot(yout(:,1), dsnum)
plot(yout(:,1), dsf, '--')
hold off
legend('dsx','dsy','dsz','fdsx','fdsy','fdsz')
