%%************************************************************************
%% soft-thresholding operator:
%% [Xp,par] = threshold(X,par,tol);
%%   Omega1 := Omega_alpha,alpha
%%   Omega2 := Omega_alpha,gamma
%%   Omega3 := Omega_alpha,beta
%%************************************************************************

function [Xp,par] = threshold2(X,par,tol)

tmpversion = version('-release');
matlabrelease = str2num(tmpversion(1:4));

if (nargin == 2); tol = 1e-12; end
[p,q] = size(X);
if isfield(par,'sig')
    rho = par.rho*par.sig;  % rho = rho*sigma;
else
    rho = par.rho;
end

if (matlabrelease < 2011)
    if (p==q)
        [U,S,V1] = mexsvd(X);
    elseif (par.use_QR)
        [U,S,V1] = mexsvd(X);
        % householder vectors v are stored in tril(HouseQ);
        % size(HouseQ) = size(V1);
        par.HouseQ = mexhouse(V1);
    else
        [U,S,V] = mexsvd(X,1);
        V1 = V(:,[1:p]);
        par.V2 = V(:,[(p+1): q]);
        par.V2t = par.V2';
    end
    %%
    d=diag(S);
    %%
    %%  The lapack function DGESDD could compute negative singular values
    %%
    z = randn(q,1);
    tmpz = U*(d.*(V1'*z));
    Xz = X*z;  normXz = norm(Xz);
    ratio = norm(tmpz - Xz)/(1+normXz);
    if isnan(ratio) || (ratio > 1e-12)
        if (p==q)
            [U,S,V1] = svd(X);
        elseif (par.use_QR)
            [U,S,V1] = svd(X,'econ'); % economic SVD
            par.HouseQ = mexhouse(V1);
        else
            [U,S,V] = svd(X);  % full SVD for matvec.m
            V1 = V(:,1:p);
            par.V2 = V(:,[(p+1): q]);
            par.V2t = (par.V2)';
        end
        if (~isnan(ratio))
            fprintf(' mexsvd ratio = %3.2e ', ratio);
        else
            fprintf(' N ');
        end
    end
else
    if (p==q)
        [U,S,V1] = svd(X);
    elseif (par.use_QR)
        [U,S,V1] = svd(X,'econ');
        % householder vectors v are stored in tril(HouseQ);
        % size(HouseQ) = size(V1);
        par.HouseQ = mexhouse(V1);
    else
        [U,S,V] = svd(X);  % full SVD for matvec.m
        V1 = V(:,[1:p]);
        par.V2 = V(:,[(p+1): q]);
        par.V2t = par.V2';
    end
end

%%
V1t = V1';
d = diag(S);
dtmp = d-rho;
posidx = find(dtmp>tol);
if isempty(posidx)
    Xp = zeros(p,q);
    dp = []; dtp = [];
    % normXp=0;    % Frobenius norm of Xp
elseif (length(posidx) == p)
    SS = spdiags(dtmp,0,p,p);
    Xp = U*SS*V1t;
    dp = d; dtp=dtmp;
    % normXp = norm(dtp);  % Frobenius norm of Xp
    %%
    Gamma = (dtmp*ones(1,p) + ones(p,1)*dtmp')./(d*ones(1,p) + ones(p,1)*d');
    par.Gamma = Gamma;
    if (p<q)
        % Beta = (dtmp*ones(1,q-p))./(d*ones(1,q-p));
        Beta = (dtmp./d)*ones(1,q-p);
        par.Beta = Beta;
    end
else
    r = length(posidx); s = p-r;
    negidx = [r+1:p];
    dp = d(posidx);
    dtp = dtmp(posidx);
    % normXp = norm(dtp);  % Frobenius norm of Xp
    dn = d(negidx);
    UU = U(:,posidx);
    VVt = V1t(posidx,:);
    SS = spdiags(dtp,0,r,r);
    Xp = UU*SS*VVt;
    %%
    tmp1= (dtp*ones(1,s))./(dp*ones(1,s) - ones(r,1)*dn');
    % Omega1 = [ones(r,r) tmp1; tmp1' zeros(s,s)];   % Omega1 = Omega_alpha,alpha
    par.Omega12 = tmp1;
    par.Omega12t = tmp1';
    tmp2= (dtp*ones(1,r) + ones(r,1)*dtp')./(dp*ones(1,r) + ones(r,1)*dp');
    tmp3= (dtp*ones(1,s))./(dp*ones(1,s) + ones(r,1)*dn');
    % Omega2=[tmp2 tmp3; tmp3' zeros(s,s)]; % Omega2 = Omega_alpha,gamma
    par.Gamma11 = tmp2;
    par.Gamma12 = tmp3;
    par.Gamma12t = tmp3';
    if (p<q)
        % tmp4 = (dtp*ones(1,q-p))./(dp*ones(1,q-p));
        tmp4 = (dtp./dp)*ones(1,q-p);
        par.Beta = tmp4;
    end
end
%%
par.U  = U;
par.V1 = V1;
par.V1t = V1t;
par.d = d;   % singular values of W
par.dp = dp;  % singular values greater than rho*sig
par.dtp = dtp; % singular values of Xp;
par.posidx = posidx;

%%************************************************************************
