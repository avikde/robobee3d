function [p,Rb,dq,nmeas] = stateEst(t, qvicon)
	% Low-pass filter
	% persistent Rbprev Rbdot eulprev v pprev tprevM
	persistent tprevP tprevM qvprev xhat Rbhat A Q R H P nmeas1 Rbvprev omgfilt lpfomg lpfR
	if isempty(tprevP)
		%fprintf(1, 'Initializing\n')
		tprevM = 0;
		tprevP = 0;
		qvprev = zeros(6,1);
		xhat = zeros(6,1); % positions
		Rbhat = eye(3);
		Rbvprev = eye(3);
		A = eye(6);
		Q = 1e-1 * diag([0.1 0.1 0.1 100 100 100]);
		R = 1e2 * diag([0.1 0.1 0.1]);
		H = [eye(3) zeros(3,3)];
		P = eye(6);
		nmeas1 = 0;
		omgfilt = zeros(3,1);
		lpfomg = 0.1;
		lpfR = 0.3;
	end
	
	% Prediction
	dt = t - tprevP;
	tprevP = t;
	A(1:3, 4:6) = dt * eye(3);
	xhat = A * xhat;
	P = A * P * A' + Q;
	% group
	Rbhat = Rbhat * expm(skew(omgfilt) * dt); % checked order; correct
	% not guaranteed SO(3) after this
	
	% Measurement if the vicon data changed
	if sum(isnan(qvicon)) == 0 && norm(qvicon - qvprev) > 1e-3
		nmeas1 = nmeas1 + 1;
		dt = t - tprevM;
		
		y = qvicon(1:3) - H * xhat;
		S = H * P * H' + R;
		K = P * H' / S;
		xhat = xhat + K * y;
		P = (eye(6) - K * H) * P;
		
		qvprev = qvicon;
		Rbv = eul2rotm(qvicon(4:6)', 'ZYZ');

		% Find Rbdot
		Rbdot = (Rbv - Rbvprev) / dt;
		Rbvprev = Rbv;
		% body-frame ang vel
		omgfilt = omgfilt + lpfomg * (unskew(Rbv' * Rbdot) - omgfilt);
		
		% Correct orientation
		errHat = -Rbv' * Rbhat + Rbhat' * Rbv;
		Rbhat = Rbhat * expm(lpfR * errHat);
		
		tprevM = t;
	end
	
	nmeas = nmeas1;
	% Outputs
	Rb = Rbhat;
	dq = [xhat(4:6); omgfilt];
	p = xhat(1:3);
end

function X = skew(v)
	X=[0   -v(3)  v(2);
	 v(3)    0  -v(1);
	-v(2)    v(1)   0  ];
end

function v = unskew(X)
	v = [X(3,2);X(1,3);X(2,1)];
end