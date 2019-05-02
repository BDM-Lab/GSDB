% rotmatfrom axis, angle
%
% c = cos(theta); s = sin(theta); C = 1-c
% xs = x*s;   ys = y*s;   zs = z*s
% xC = x*C;   yC = y*C;   zC = z*C
% xyC = x*yC; yzC = y*zC; zxC = z*xC
% [ x*xC+c   xyC-zs   zxC+ys ]
% [ xyC+zs   y*yC+c   yzC-xs ]
% [ zxC-ys   yzC+xs   z*zC+c ]
%
% http://en.wikipedia.org/wiki/Rotation_matrix
%
% Q = rotmatfrom(myaxis,myangle)

function Q = rotmatfrom(myaxis,myangle)

x = myaxis(1);
y = myaxis(2);
z = myaxis(3);
theta = myangle;

c = cos(theta);
s = sin(theta);
C = 1-c;
xs = x*s;   ys = y*s;   zs = z*s;
xC = x*C;   yC = y*C;   zC = z*C;
xyC = x*yC; yzC = y*zC; zxC = z*xC;

Q = [ x*xC+c   xyC-zs   zxC+ys ;
      xyC+zs   y*yC+c   yzC-xs ;
      zxC-ys   yzC+xs   z*zC+c];

