% This is all a guess
function [h0, pdotdes] = mpc2wl(accdes, Rb, mb, Ib, g)
	h0 = [Rb' * [0; 0; mb*g]; 0;0;0]
	pdotdes = diag([mb,mb,mb,Ib]) * accdes
end

% ---------
function actualT0 = wl2mpc(w0, mb)
	actualT0 = w0[2]/mb
end
