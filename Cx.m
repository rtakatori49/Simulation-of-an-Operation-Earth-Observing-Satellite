function R = Cx(theta)
% R1  Find the rotation matrix about the number one axis for a rotation of
% theta radians
%

R = [1 0          0;
     0 cos(theta) sin(theta);
     0 -sin(theta) cos(theta)];