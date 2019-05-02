%%************************************************************************
%% soft-thresholding operator:
%% [Xp,par] = threshold(X,par,tol);
%%   Omega1 := Omega_alpha,alpha
%%   Omega2 := Omega_alpha,gamma
%%   Omega3 := Omega_alpha,beta
%%************************************************************************

    function [Xp,par] = threshold2old(X,par,tol)
    
    if (nargin == 2);
       tol = 1e-12;
    end
    [p,q] = size(X);
    if isfield(par,'sig')
        rho = par.rho*par.sig;  % rho = rho*sigma;
    else
        rho = par.rho;
    end
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
    V1t = V1';
    d=diag(S);
    dtmp= d-rho;
    posidx = find(dtmp>tol);
    if isempty(posidx)
       Xp = zeros(p,q);
       Omega1 = zeros(p,p); Omega2 = zeros(p,p); 
       if (p<q)
           Omega3 = zeros(p,q-p);
           par.Omega3 = Omega3;
       end
       dp=[]; dtp=[];
       % normXp=0;    % Frobenius norm of Xp
    elseif (length(posidx) == p)
       SS = spdiags(dtmp,0,p,p);
       Xp = U*SS*V1t;
       Omega1 = ones(p,p);
       Omega2 = (dtmp*ones(1,p) + ones(p,1)*dtmp')./(d*ones(1,p) + ones(p,1)*d');
       if (p<q)
           % Omega3 = (dtmp*ones(1,q-p))./(d*ones(1,q-p));
           Omega3 = (dtmp./d)*ones(1,q-p);
           par.Omega3 = Omega3;
       end
       dp = d; dtp=dtmp;
       % normXp = norm(dtp);  % Frobenius norm of Xp
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
       Omega1 = [ones(r,r) tmp1; tmp1' zeros(s,s)];   % Omega1 = Omega_alpha,alpha
       tmp2= (dtp*ones(1,r) + ones(r,1)*dtp')./(dp*ones(1,r) + ones(r,1)*dp');
       tmp3= (dtp*ones(1,s))./(dp*ones(1,s) + ones(r,1)*dn');
       Omega2=[tmp2 tmp3; tmp3' zeros(s,s)]; % Omega2 = Omega_alpha,gamma
       if (p<q)
           % tmp4 = (dtp*ones(1,q-p))./(dp*ones(1,q-p));
           tmp4 = (dtp./dp)*ones(1,q-p);
           Omega3 = [tmp4; zeros(s,q-p)];   % Omega3 = Omega_alpha,beta
           par.Omega3 = Omega3;
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
    par.Omega1=Omega1;
    par.Omega2=Omega2;
    % par.Omega3=Omega3;

%%************************************************************************
