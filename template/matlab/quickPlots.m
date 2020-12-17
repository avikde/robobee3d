function quickPlots(yout, sampling_time, tmax)

clf
subplot(321)
plot(yout(:,1), yout(:,13))
xlim([0,tmax])
ylabel('V')
% 
% if 0
subplot(322)
hold all
plot(yout(:,1), yout(:,14))
plot(yout(:,1), yout(:,15))
legend('drv_pch', 'drv_roll')
xlim([0,tmax])
% end 
subplot(323)
hold all
plot(yout(:,1), qvicon(:,1:3))
plot(yout(:,1), plog, '--')
if size(yout,2) > 34
	posdes = yout(:,35:37);
	plot(yout(:,1), posdes(:,2), 'k--')
	legend('raw x','raw y','raw z', 'f x','f y','f z', 'desy')
else
	legend('raw x','raw y','raw z', 'f x','f y','f z')
end
ylabel('pos')
xlim([0,tmax])


subplot(324)
% hold all
% plot(yout(:,1), qvicon(:,4:6))
% ylim([-1,1])
% ylabel('rot vicon')
% legend('r rx','r ry','r rz')

hold all
Np=size(yout,1);
omgnum=zeros(Np,3);
for i=2:Np
	o0 = omgnum(i-1,:);
	o1 = (qvicon(i,4:6) - qvicon(i-1,4:6))/sampling_time;
	omgnum(i,:) = o0 + 0.005*(o1 - o0);
end
plot(yout(:,1), dqlog(:,4:6))
%plot(yout(:,1), omgnum, '--')
ylabel('omega')
legend('wx','wy','wz')%,'nwx','nwy','nwz')

Np = size(yout,1);
ss = zeros(Np, 3);
dss = zeros(Np, 3);
e3h = [0 -1 0; 1 0 0; 0 0 0];
% approximate if caught on tether
panchor = [0,-0.05,25*0.0254];
ltether = zeros(Np,1);
for i = 1:Np
	Rb = eul2rotm(qvicon(i,4:6),'ZYZ');
% 	Rb = get_rotation_matrix(qvicon(i,4:6));
	ss(i,:) = Rb * [0; 0; 1];
	dss(i,:) = (-Rb * e3h * dqlog(i,4:6)')';
	
	ltether(i) = norm(plog(i,:) - panchor);
end
ltether = (1 - ltether/ltether(1)) * 5;
xlim([0,tmax])

subplot(325)
% hold all
% plot(yout(:,1), ss)
% %plot(yout(:,1), ltether, 'k')
% hold off
% legend('sx','sy','sz')%,'lt')
% 
hold all
plot(yout(:,1), uquadlog(:,1), 'r--')
plot(yout(:,1), 1e-3*uquadlog(:,2), 'g')
plot(yout(:,1), 1e-3*uquadlog(:,3), 'b')
legend('thrust','momx','momy')
%ylabel('V')
hold off
xlim([0,tmax])

subplot(326)
hold all
plot(yout(:,1), dss)
plot([yout(1,1) yout(end,1)], [0 0], 'k--')
ylim([-30,30])
hold off
legend('dsx','dsy','dsz')
xlim([0,tmax])
