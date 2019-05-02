% rotmatto axis, angle
%
% x = Qzy-Qyz
% y = Qxz-Qzx
% z = Qyx-Qxy
% r = sqrt(x^2+y^2+z^2)
% t = Qxx+Qyy+Qzz
% theta = atan2(r,t-1)
%
% Reference: http://en.wikipedia.org/wiki/Rotation_matrix
%
% [myaxis,myangle] = rotmatto(Q)

function [myaxis,myangle] = rotmatto(Q)

x = Q(3,2)-Q(2,3);
y = Q(1,3)-Q(3,1);
z = Q(2,1)-Q(1,2);
r = sqrt(x^2+y^2+z^2);
t = trace(Q);
theta = atan2(r,t-1);

myaxis = [x y z]' / r;
myangle = theta;

