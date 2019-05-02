% computermsd Compute the RMSD between two sets of points.
%
% rmsd = computermsd(A,B,shouldAlign)

% Created       :  6 Jul 07
% Last modified :  5 Nov 07

function rmsd = computermsd(A,Aest,shouldAlign)

if shouldAlign
    Info = alignatoms(A,Aest);
    Aest = Info.Aest;
end

n    = size(A,2);
minx=min(min(A));%min(min(min(A)),min(min(Aest)));
maxx=max(max(A));%max(max(max(A)),max(max(Aest)));
rmsd = norm(A-Aest,'fro') / sqrt(n)/(maxx-minx);

