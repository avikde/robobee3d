function wltest()

% repeatedly call for a made-up pdes to test
Nn=500;
aa = zeros(Nn,3);
for i = 1:Nn
	[drv_amp, drv_pch, drv_roll] = wlwrap(zeros(6), [0,0,10,0,0,0]);
	aa(i,:) = [drv_amp, drv_pch, drv_roll];
end

figure(1)
subplot(2,1,1)
plot(aa(:,1))
subplot(2,1,2)
plot(aa(:,2:3))

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
pddes = [0, 0, 10, 0,0,0];

% Update
u = wlControllerUpdate(single(u0), single(h0), single(pddes));

% Convert to robot coords. u = [Vmean, uoffs, udiff, h2]. uoffs, udiff
% normalized by Vmean
drv_amp = u(1);
drv_pch = drv_amp * u(2);
drv_roll = drv_amp * u(3);

u0 = u;

end
