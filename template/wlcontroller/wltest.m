function wltest()

% repeatedly call for a made-up pdes to test
Nt = 500;
U = zeros(Nt,4);
pdotdes = [0,0,10,0,0,0];

U(1,:) = [140.0, 0., 0., 0.];
mb = 100;
g = 9.81e-3;
h0 = [0, 0, mb * g, 0, 0, 0];

for i = 2:Nt
	if i > 100
			pdotdes(3) = 20;
	end
	if i > 200
			pdotdes(3) = 30;
	end
	if i > 300
			pdotdes(3) = 40;
	end
	if i > 400
			pdotdes(3) = 50;
	end
	U(i,:) = wlControllerUpdate(single(U(i-1,:)), single(h0), single(pdotdes));
end

% convert
Vmean = U(:,1);
uoffs = U(:,2);
udiff = U(:,3);

vleft = Vmean .* (1 + udiff) .* (1 + uoffs);
vright = Vmean .* (1 - udiff) .* (1 + uoffs);

drv_pch = Vmean .* uoffs;
drv_amp = abs((vleft-vright))/2+min(vleft,vright);
drv_roll = (vleft-vright)/4; 
drv_bias = max(vleft,vright)+2*abs(drv_pch); % voltage

subplot(2,1,1)
hold all
plot(Vmean)
plot(drv_bias)
plot(drv_amp)
ylim([100,220])
legend('Vmean','drv_bias','drv_amp')

subplot(2,1,2)
hold all
plot(drv_roll)
plot(drv_pch)
ylim([-30,30])
legend('drv_roll','drv_pch')

