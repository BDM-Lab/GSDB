% isomatto sign, axis, angle
%
% [mysign,myaxis,myangle] = isomat2(Q)

function [mysign,myaxis,myangle] = isomatto(Q)

mysign = sign(det(Q));
[myaxis,myangle] = rotmatto(mysign*Q);

