% rotate Rotate points by a rotation matrix.
%
% points = rotate(points,rotationMatrix)

% Created       :  6 Jul 07
% Last modified :  6 Jul 07

function points = rotate(points,rotationMatrix)

points = rotationMatrix * points;
