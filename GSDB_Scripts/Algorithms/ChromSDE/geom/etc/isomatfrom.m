% isomatfrom sign, axis, angle
%
% Q = isomatfrom(mysign,myaxis,myangle)

function Q = isomatfrom(mysign,myaxis,myangle)

Q = mysign * rotmatfrom(myaxis,myangle);

