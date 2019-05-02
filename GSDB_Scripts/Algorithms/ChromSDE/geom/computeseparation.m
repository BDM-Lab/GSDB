% computeseparation Compute the separation of a set of points.
%
% The separation is defined as
%                          sqrt(sum(||x_i-x_j||^2)).
%
% separation = computeseparation(points)

% Created       :  6 Jul 07
% Last modified :  6 Jul 07

function separation = computeseparation(points)

nPoints    = size(points,2);
separation = 0;
for iPoint = 1:nPoints
    pointsij   = points(:,1:(iPoint-1)) - points(:,iPoint) * ones(1,iPoint-1);
    separation = separation + sum(sum(pointsij.^2));
end
separation = sqrt(separation);

