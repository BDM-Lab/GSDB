%%************************************************************************
%% project2: compute projection onto the cone of positive
%%           semidefinite matrices.
%%
%% [Xp,par] = project2(blk,X,par);
%%
%% SDPNAL:
%% Copyright (c) 2008 by
%% Xinyuan Zhao, Defeng Sun, and Kim-Chuan Toh
%%************************************************************************

function [Xp,par] = project2(blk,X,par);

tol = 1e-12;  %% must be small as it will affect gap
addtol = 1e-6;
Xp = cell(size(X));
par.P1 = cell(size(X));     par.P2 = cell(size(X));
par.Dsch12 = cell(size(X)); par.dd = cell(size(X));
par.posidx = cell(size(X));
%%
for p = 1:size(blk,1)
    pblk = blk(p,:);
    if strcmp(pblk{1},'s');
        blktmp = pblk{1,2};
        if (length(blktmp) > 1);
            %%error(' each cell can only have one block');
        end
        n = sum(blktmp);
        [Xptmp,P,Dsch12,dd,posidx] = projectfun2(X{p},tol);
        %%***** perturbation *****
        Dsch12 = max(addtol,Dsch12);
        if ~isempty(posidx)
            P1 = P(:,posidx);
            P2 = P(:,setdiff([1:n],posidx));
        else
            P1 = [];
            P2 = P;
        end
    elseif strcmp(pblk{1},'l');
        P1 = []; P2 = [];
        n  = sum(pblk{1,2});
        Xptmp = zeros(n,1);
        %%***** perturbation *****
        Dsch12 = addtol*ones(n,1);
        posidx = find(X{p} > tol);
        if ~isempty(posidx)
            Xptmp(posidx)  = abs(X{p}(posidx));
            Dsch12(posidx) = ones(length(posidx),1);
        end
        dd = [];
    end
    Xp{p}     = Xptmp;
    par.P1{p} = P1;
    par.P2{p} = P2;
    par.dd{p} = dd;
    par.posidx{p} = posidx;
    par.Dsch12{p} = Dsch12;
end
%%***************************************************************************

function [Xp,V,Dsch12,d,posidx] = projectfun2(X,tol);

n = length(X);
if (norm(X-X','fro') > 1e-15*norm(X,'fro'))
    fprintf(' nonSym ');
    X = 0.5*(X+X');
end
if (exist('mexeig')==3)
    [V,D] = mexeig(full(X));
else
    [V,D] = eig(full(X));
end
if (size(D,1)==size(D,2))
   d = diag(D);
end
[d,idx] = sort(real(d));
idx = idx(n:-1:1); d = d(n:-1:1);
V = V(:,idx);
posidx = find(d > tol);
if isempty(posidx)
    Xp = sparse(n,n);
    Dsch12 = [];
elseif (length(posidx) == n)
    Xp = X;
    Dsch12 = [];
else
    r = length(posidx); s = n-r;
    negidx = [r+1:n];
    dp = abs(d(posidx));
    dn = abs(d(negidx));
    Vtmp = V(:,posidx)*diag(sqrt(dp));
    Xp = Vtmp*Vtmp';
    Xp = 0.5*(Xp+Xp');
    Dsch12 = (dp*ones(1,s))./(dp*ones(1,s) + ones(r,1)*dn');
end
%%************************************************************************
