function wltest()

% repeatedly call for a made-up pdes to test

end

% ---------------------

function [drv_amp, drv_pch, drv_roll] = wlwrap(viconState, pdes)

% Store u0 to update
persistent u0
if isempty(u0)
	u0 = zeros(4);
end


u0 = [140.0, 0., 0., 0.];
mb = 100;
g = 9.81e-3;
h0 = [0, 0, mb * g, 0, 0, 0];
p0 = zeros(1,6);
pdes = [0, 0, 20, 0, 0, 0];
Qdiag = [1,1,1,0.1,0.1,0.1];
kpmom = [0,0,1,0.1,0.1,0.1];

u = wlmex(u0, p0, h0, pdes, kpmom, Qdiag);

drv_amp = 150;
drv_pch = 0;
drv_roll = 0;

end
