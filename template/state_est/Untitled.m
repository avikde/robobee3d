load('../../../../Desktop/20200926_closedloop/p10_d0_2.mat')
set(0, 'DefaultLineLineWidth', 2)
Np = size(yout, 1);

% 
s = zeros(Np,3);
ds = zeros(Np,3);
dsnum = zeros(Np,3);
omgnum = zeros(Np,3); % from Rb
ww = zeros(Np, 3); % differentiate ZYZ angles NOT WORKING
e3h = [0 -1 0; 1 0 0; 0 0 0];
Rbprev = eye(3);
for i=1:Np
	Rb = eul2rotm(yout(i,10:12), 'ZYZ');
	s(i,:) = Rb * [0 0 1]';
	if i > 1
		dt = yout(i,1) - yout(i-1,1);
		sdiff = (s(i,:) - s(i-1,:)) / dt;
		Rbdot = (Rb - Rbprev) / dt;
		omghat = Rb' * Rbdot;
		omg1 = [omghat(3,2),omghat(1,3),omghat(2,1)];
		omgnum(i,:) = omgnum(i-1,:) + 0.005 * (omg1 - omgnum(i-1,:));
		dsnum(i,:) = dsnum(i-1,:) + 0.005 * (sdiff - dsnum(i-1,:));
		%
		angdiff = (yout(i,19:21) - yout(i-1,19:21)) / dt;
		ww(i,:) = ww(i-1,:) + 0.005 * (angdiff - ww(i-1,:));
	end
	ds(i,:) = -Rb * e3h * omgnum(i,:)';
	Rbprev = Rb;
end

clf
subplot(221)
hold all
plot(yout(:,1), yout(:,10:12))
plot(yout(:,1), yout(:,19:21), '--')
hold off
legend('rx','ry','rz','frx','fry','frz')

subplot(222)
hold all
plot(yout(:,1), omgnum, '--')
plot(yout(:,1), yout(:,25:27))
hold off
legend('drx','dry','drz','fdrx','fdry','fdrz')


subplot(223)
hold all
%plot(yout(:,1), s)
%legend('sx','sy','sz')
plot(yout(:,1), ww, '--')
plot(yout(:,1), omgnum)
hold off
legend('drx','dry','drz','ndrx','ndry','ndrz')

subplot(224)
hold all
plot(yout(:,1), ds)
plot(yout(:,1), dsnum, '--')
hold off
legend('dsx','dsy','dsz','ndsx','ndsy','ndsz')
