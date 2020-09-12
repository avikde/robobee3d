
load('../../../../Dropbox (Personal)/harvard/micro/simulink/wlcontroller/teststateest.mat')
size(yout)

%%

Nt = size(yout,1);
qfilt = zeros(Nt, 6);
dqfilt = zeros(Nt, 6);

for i=1:Nt
	[q, Rb, dq] = stateEst(yout(i,7:12)',10000,500);
	qfilt(i,:) = q;
	dqfilt(i,:) = dq;
end

subplot(221)
hold all
plot(yout(:,1), yout(:,7:9), '.')
plot(yout(:,1), qfilt(:,1:3))
ylabel('pos')
legend('x','y','z')

subplot(222)
hold all
plot(yout(:,1), yout(:,10:12), '.')
plot(yout(:,1), qfilt(:,3:5))
ylabel('rot')
legend('rx','ry','rz')

subplot(223)
hold all
plot(yout(:,1), dqfilt(:,1:3))
ylabel('vel')
legend('x','y','z')

subplot(224)
hold all
plot(yout(:,1), dqfilt(:,3:5))
ylabel('velrot')
legend('rx','ry','rz')
