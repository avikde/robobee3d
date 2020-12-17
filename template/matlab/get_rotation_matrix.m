function [R]  = get_rotation_matrix(theta)

alpha = theta(1); beta = theta(2); gamma = theta(3);
x = alpha; y = beta; z = gamma;

R11 = cos(y)*cos(z);
R12 = -cos(y)*sin(z);
R13 = sin(y);

R21 = cos(x)*sin(z)+sin(x)*sin(y)*cos(z);
R22 = cos(x)*cos(z)-sin(x)*sin(y)*sin(z);
R23 = -sin(x)*cos(y);

R31 = sin(x)*sin(z)-cos(x)*sin(y)*cos(z);
R32 = sin(x)*cos(z)+cos(x)*sin(y)*sin(z);
R33 = cos(x)*cos(y);

R = [R11 R12 R13; R21 R22 R23; R31 R32 R33];

end