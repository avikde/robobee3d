function wltest()

% repeatedly call for a made-up pdes to test
aa = zeros(100,3);
for i = 1:100
	[drv_amp, drv_pch, drv_roll] = wlwrap(zeros(6), [0,0,10,0,0,0]);
	aa(i,:) = [drv_amp, drv_pch, drv_roll];
end

figure(1)
for i = 1:3
	subplot(3,1,i)
	plot(aa(:,i))
end

end

% ---------------------

function [drv_amp, drv_pch, drv_roll] = wlwrap(viconState, pdes)

% Store u0 to update
persistent u0
if isempty(u0)
	u0 = [140.0, 0., 0., 0.];
end

% u0 = [140.0, 0., 0., 0.];
mb = 100;
g = 9.81e-3;
h0 = [0, 0, mb * g, 0, 0, 0];
p0 = zeros(1,6);
% pdes = [0, 0, 20, 0, 0, 0];
Qdiag = [1,1,1,0.1,0.1,0.1];
kpmom = [0,0,1,0.1,0.1,0.1];

% Update
u = wlmex(u0, p0, h0, pdes, kpmom, Qdiag);

% Convert to robot coords. u = [Vmean, uoffs, udiff, h2]. uoffs, udiff
% normalized by Vmean
drv_amp = u(1);
drv_pch = drv_amp * u(2);
drv_roll = drv_amp * u(3);

u0 = u;

end
