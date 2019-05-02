
%% diagonal preconditioner for IDR(s)
%% L.invdiagM

function L = diagapprox(AAt,par)

p = par.p;  q = par.q;
r = length(par.posidx);
if (r==0)
    diagH = max(1e-4,par.msk);
    L.invdiagM = 1./diagH;
    return
end

Usquare = par.U.*par.U;
% Vsquare = par.V1.*par.V1;
V1tsquare = par.V1t.*par.V1t;
% M = 0.5*(par.Omega1+par.Omega2);
if (p<q)
    % M = [M  par.Omega3];
    if (par.use_QR)
        tmp = eye(q-p,q-p);
        V2t = mexAV2(tmp,par.HouseQ,1);
        V2tsquare = V2t.*V2t;
        Vtsquare = [V1tsquare; V2tsquare];
    else
        V2tsquare = par.V2t.*par.V2t;
        Vtsquare = [V1tsquare; V2tsquare];
    end
else
    Vtsquare = V1tsquare;
end

if (r==p)
    M = 0.5*(ones(p,p)+par.Gamma);
    if (p<q)
        M = [M par.Beta];
    end
    D = M*Vtsquare;
    D = Usquare*D;
else
    M11 = 0.5*(ones(r,r)+par.Gamma11);
    M12 = 0.5*(par.Omega12 + par.Gamma12);
    M21 = 0.5*(par.Omega12t + par.Gamma12t);
    tmpM = [M11  M12];
    if (p<q)
        M13 = par.Beta;
        tmpM = [tmpM  M13];
    end
    tmp1 = tmpM*Vtsquare;
    % V1rt = par.V1t([1:r],:);
    % V1rtsquare = V1rt.*V1rt;
    V1rtsquare = V1tsquare([1:r],:);
    tmp2 = M21*V1rtsquare;
    tmp3 = [tmp1; tmp2];
    D = Usquare*tmp3;
end
%%
d = D(:);
diagAAt = (d'*AAt)';  % AAt = At.*At;
diagAAt = max(1e-4,full(diagAAt));
if isfield(par,'sig'); diagAAt = par.sig*diagAAt; end;
diagH = par.msk+ diagAAt;
L.invdiagM = 1./diagH;


