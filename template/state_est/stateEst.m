function [q,Rb,dq,nmeas] = stateEst(t, qvicon)
	% Low-pass filter
	% persistent Rbprev Rbdot eulprev v pprev tprevM
	persistent tprevP tprevM qvprev xhat Rbhat A Q R H P nmeas1 Rbvprev omgfilt
	if isempty(tprevP)
		%fprintf(1, 'Initializing\n')
% 		Rbprev = eye(3);
% 		Rbdot = zeros(3,3);
% 		eulprev = zeros(3,1);
% 		pprev = zeros(3,1);
% 		v = zeros(3,1);
		tprevM = 0;
		tprevP = 0;
		qvprev = zeros(6,1);
		xhat = zeros(12,1);
		Rbhat = eye(3);
		Rbvprev = eye(3);
		A = eye(12);
		Q = 1e-1 * diag([0.1 0.1 0.1 1 1 1 100 100 100 1000 1000 1000]);
		R = 1e2 * diag([0.1 0.1 0.1 1 1 1]);
		H = [eye(6) zeros(6,6)];
		P = eye(12);
		nmeas1 = 0;
		omgfilt = zeros(3,1);
	end
	
	% Prediction
	dt = t - tprevP;
	tprevP = t;
	%size(A(7:12, 1:6))
	A(1:6, 7:12) = dt * eye(6);
	xhat = A * xhat;
	P = A * P * A' + Q;
	% group
	%Rbhat = Rbhat + (Rbhat * skew(dqhat(3:5))) * dt;
	% not guaranteed SO(3) after this
	
	% Measurement if the vicon data changed
	if sum(isnan(qvicon)) == 0 && norm(qvicon - qvprev) > 1e-3
		nmeas1 = nmeas1 + 1;
		dt = t - tprevM;
		
		y = qvicon - H * xhat;
		S = H * P * H' + R;
		K = P * H' / S;
		xhat = xhat + K * y;
		P = (eye(12) - K * H) * P;
		
		qvprev = qvicon;
		Rbv = eul2rotm(qvicon(4:6)', 'ZYZ');

		% Find Rbdot
		Rbdot = (Rbv - Rbvprev) / dt;
		Rbvprev = Rbv;
		% body-frame ang vel
		omgHat = Rbv' * Rbdot;
		omgfilt = omgfilt + 0.1 * ([omgHat(3,2);omgHat(1,3);omgHat(2,1)] - omgfilt);
% 		
% 		eulprev = eulprev + lpfomega * ([omgHat(3,2);omgHat(1,3);omgHat(2,1)] - eulprev);
% 		% vel
% 		v = v + lpfv * ((qvicon(1:3) - pprev) * dt - v);
% 		pprev = qvicon(1:3);
		
		tprevM = t;
	end
	
	dq = xhat(7:12);
	q = xhat(1:6);
	Rb = eul2rotm(q(4:6)', 'ZYZ');
	nmeas = nmeas1;
	%TEST
	xhat(10:12) = omgfilt;
end

% function X = skew(v)
% 	X=[0   -v(3)  v(2);
% 	 v(3)    0  -v(1);
% 	-v(2)    v(1)   0  ];
% end