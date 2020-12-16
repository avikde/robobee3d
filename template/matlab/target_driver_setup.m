clear
sampling_f = 5000;			% controller sampling rate
sampling_time = 1/sampling_f;
start_delay = 0.0;		%2.0	% start after xx s

%% Controller setup

pos0 = [0.05,-0.02,0.08];
trajFreq = 0.5;
trajAmp = 0.07;

kquad = [4e1, 1e-4, 4e-5, 153, 1e0]; %kPz, kmomy, kmomx, V0, kDz
% ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom
wmpc = single([3e2, 5e1, 1e-4, 2e-4, 1e0, 1e0, 1e-1, 1e0]);
lpfParams = [0.2, 0.5, 1e-1, 1e1]; %lpfOmg, lpfR, kfQ, kfR. was 0.2, 0.5, 1e-1, 1e2

% WL inv dyn controller (mod3 file)
popts = single([-5.22e-01, 7.45e-03, 9.72e-01,-3.99e-02, 1.00e+00,-5.27e-05,-6.76e-03, 7.85e-04, 1.00e+00,-5.22e-01, 1.77e+00, 1.00e+00,-1.40e+00, 1.00e+00, 1.00e+00, 3.12e-03,-1.48e-04,-8.33e-02, 1.46e-02, 1.00e+00, 1.60e-06, 4.51e-04,-1.24e-04, 1.00e+00, 3.83e-01,-1.24e-01, 1.00e+00, 5.65e-02, 1.00e+00, 1.00e+00, 3.23e+01,-4.40e-01,-2.02e+00, 6.84e+01, 1.00e+00, 3.08e-03, 1.55e-02,-4.66e-01, 1.00e+00,-1.00e+01,-1.73e+00, 1.00e+00, 7.78e+01, 1.00e+00, 1.00e+00, 2.26e+02,-3.08e+00,-1.25e+01, 6.23e+02, 1.00e+00, 2.13e-02, 8.82e-02,-4.37e+00, 1.00e+00,-4.62e+01,-1.53e+01, 1.00e+00, 5.40e+02, 1.00e+00, 1.00e+00, 5.66e+00,-7.64e-02, 3.19e+01, 2.44e+01, 1.00e+00, 5.24e-04,-2.42e-01,-1.55e-01, 1.00e+00,-5.95e+00, 2.32e+01, 1.00e+00, 1.30e+01, 1.00e+00, 1.00e+00, 3.51e-01,-7.50e-03,-2.13e+01, 1.67e+01, 1.00e+00, 7.47e-05, 1.53e-01,-1.15e-01, 1.00e+00,-3.11e+01,-1.27e+01, 1.00e+00, 1.66e+00, 1.00e+00, 1.00e+00]);

%% setup code 
f = 160; %165
f_int = f;					% initial frequency
f_fin = f;					% final frequency
T = 2;			% total time in seconds

if 0    % for zhi
    drv_bias = 200;
    zamp = 180;

    drv_amp = zamp/2;
    drv_pch = 0;
    drv_roll = zamp/4;

else
    %v = 220;%180;
    vleft = 170; %205;%+13-5; %160 %205
    vright = 170; %193;%-13+5;  %200
    drv_pch = 0; %15;%21.5+30;%22;
    drv_amp = abs((vleft-vright))/2+min(vleft,vright);
    drv_roll = (vleft-vright)/4; 
    drv_bias = max(vleft,vright)+2*abs(drv_pch); % voltage
    phase = 0.0; %degrees (added on left wing)

%     % to rail (+) -- actuator testing
%     drv_amp = v;
%     drv_bias = v;
%     drv_pch = v;
    
%     % to rail (0) -- actuator testing
%     drv_amp = 0;
%     drv_bias = v;
%     drv_pch = -v;
end

%% end of setup code
running_time = 0.0+T+0.0;   %3.0+T+0.5 %11
