function [s, ds] = uprightCalcs(eul, omega)

Np = size(eul,1);
s = zeros(Np, 3);
ds = zeros(Np, 3);
e3h = [0 -1 0; 1 0 0; 0 0 0];

for i = 1:Np
% 	Rb = eul2rotm(qvicon(i,4:6),'ZYZ');
	Rb = get_rotation_matrix(eul(i,:));
	s(i,:) = Rb * [0; 0; 1];
	ds(i,:) = (-Rb * e3h * omega(i,:)')';
end

end