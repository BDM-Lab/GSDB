%%***************************************************************************
%% project: compute X = Xp-Xn, where Xp,Xn are positive semidefinite.
%%
%% [Xp,Xn,d,rankXp] = project(blk,X);
%%
%% SDPNAL:
%% Copyright (c) 2008 by
%% Xinyuan Zhao, Defeng Sun, and Kim-Chuan Toh
%%***************************************************************************

function [Xp,Xn,d,rankXp] = project(blk,X);

tol = 1e-12; %% must be small as it will affect gap
if ~iscell(X);
    if strcmp(blk{1},'s');
        blktmp = blk{1,2};
        n = sum(blktmp);
        if (length(blktmp) > 1);
            Xp = sparse(n,n);
            Xn = sparse(n,n);
            d  = zeros(n,1);
            rankXp = 0;
            s  = [0, cumsum(blktmp)];
            for i = 1:length(blktmp)
                pos = [s(i)+1 : s(i+1)];
                [Xpsub,Xnsub,dsub,rankXpsub] = projectfun(full(X(pos,pos)),tol);
                d(pos)  = dsub;
                Xp(pos,pos) = sparse(Xpsub);
                Xn(pos,pos) = sparse(Xnsub);
                rankXp = rankXp + rankXpsub;
            end
        else
            [Xp,Xn,d,rankXp] = projectfun(X,tol);
        end
    elseif strcmp(blk{1},'l');
        n = sum(blk{1,2});
        d = sparse(n,1);
        Xp = zeros(n,1); posidx = find(X > tol);
        if ~isempty(posidx)
            Xp(posidx) = abs(X(posidx)); Xn = abs(Xp-X);
        else
            Xn = abs(X);
        end
        rankXp = length(posidx);
    end
else
    Xp = cell(size(X)); Xn = cell(size(X)); d = cell(size(X));
    for p = 1:size(blk,1)
        [Xp{p},Xn{p},d{p},rankXp(p)] = project(blk(p,:),X{p});
    end
end
%%***************************************************************************

function [Xp,Xn,d,rankXp] = projectfun(X,tol);

n = length(X);
if (norm(X-X','fro') > 1e-15*norm(X,'fro'))
    fprintf('nonSym');
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
posidx = find(d > tol);
rankXp = length(posidx);
if isempty(posidx)
    Xp = sparse(n,n); Xn = -X;
elseif (length(posidx) == n)
    Xp = X; Xn = sparse(n,n);
else
    if (length(posidx) <= n/2)
        dp = abs(d(posidx));
        Vtmp = V(:,posidx)*diag(sqrt(dp));
        Xp = Vtmp*Vtmp'; Xp = 0.5*(Xp+Xp');
        Xn = Xp-X;
    else
        negidx = setdiff([1:n],posidx);
        dn = abs(d(negidx));
        Vtmp = V(:,negidx)*diag(sqrt(dn));
        Xn = Vtmp*Vtmp'; Xn = 0.5*(Xn+Xn');
        Xp = X+Xn;
    end
end
%%***************************************************************************
