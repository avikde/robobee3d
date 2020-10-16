
%load('../../../../Dropbox (Personal)/harvard/micro/simulink/wlcontroller/teststateest.mat')
load('../../../../Desktop/20200926_closedloop/p10_d0_2.mat')
size(yout)

%%

clear stateEst

Nt = size(yout,1);
qfilt = zeros(Nt, 6);
dqfilt = zeros(Nt, 6);

for i=1:Nt
	[q, Rb, dq, nmeas] = stateEst(yout(i,1), yout(i,7:12)');
	qfilt(i,:) = q;
	dqfilt(i,:) = dq;
end
nmeas/Nt * 10000 % approx vicon rate discarding nan

subplot(221)
hold all
plot(yout(:,1), yout(:,7:9), '.')
plot(yout(:,1), qfilt(:,1:3), 'linewidth', 2)
ylabel('pos')
legend('x','y','z', 'xf','yf','zf')

subplot(222)
hold all
plot(yout(:,1), yout(:,10:12), '.')
plot(yout(:,1), qfilt(:,4:6), 'linewidth', 2)
ylabel('rot')
legend('rx','ry','rz','rxf','ryf','rzf')

subplot(223)
hold all
plot(yout(:,1), dqfilt(:,1:3), 'linewidth', 2)
ylabel('vel')
legend('x','y','z')

subplot(224)
hold all
plot(yout(:,1), dqfilt(:,4:6), 'linewidth', 2)
ylabel('velrot')
legend('rx','ry','rz')
