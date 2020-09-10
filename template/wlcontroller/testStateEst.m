function [q,Rb,dq] = stateEst(qvicon, controlRate)

[q,Rb,dq] = lowPassFilter(qvicon, controlRate);

end

function [q,Rb,dq] = lowPassFilter(qvicon, controlRate)
	% Low-pass filter
	persistent Rbprev Rbdot omega v pprev
	if isempty(Rbprev)
			Rbprev = eye(3);
			Rbdot = zeros(3,3);
			omega = zeros(1,3);
			pprev = zeros(1,3);
			v = zeros(1,3);
	end
	
	% LPF params 0 to 1. 1 => deadbeat
	lpfR = 0.1;
	lpfomega = 0.1;
	lpfv = 0.1;

	eul = qvicon(3:5);
	Rb = eul2rotm(eul');

	% Find Rbdot
	Rbdot = Rbdot + lpfR * ((Rb - Rbprev) * controlRate - Rbdot);
	Rbprev = Rb;
	% body-frame ang vel
	omgHat = Rb' * Rdot;
	omega = omega + lpfomega * ([omgHat(3,2), omgHat(1,3), omgHat(2,1)] - omega);
	% vel
	v = v + lpfv * ((qvicon(1:3) - pprev) * controlRate - v);
	pprev = qvicon(1:3);
	
	dq = [v,omega];
	q = qvicon;
end
