
load('teststateest.mat')

%%

subplot(221)
hold all
plot(yout(:,1), yout(:,7:9), '.')
ylabel('pos')
legend('x','y','z')

subplot(222)
hold all
plot(yout(:,1), yout(:,10:12), '.')
ylabel('rot')
legend('rx','ry','rz')
