%%*****************************************************************
%% diagapprox: compute an approximation to
%%             diag (A diag((PxP) Dsch (PtxPt)) At)
%%
%% SDPNAL:
%% Copyright (c) 2008 by
%% Xinyuan Zhao, Defeng Sun, and Kim-Chuan Toh
%%*****************************************************************

function dd = diagapproxsym(blk,At,par)

m = size(At{1},2);
dd = zeros(m,1);
for p = 1:size(blk,1)
    pblk = blk(p,:);
    n = sum(pblk{2});
    if strcmp(pblk{1},'s')
        PP1 = (par.P1{p}).*(par.P1{p});
        PP2 = (par.P2{p}).*(par.P2{p});
        rr = size(par.P1{p},2);
        if (rr == n)
            sumPP1 = PP1*ones(rr,1);
            tmp = sumPP1*sumPP1';
        elseif (rr > 0)
            tmp2 = (PP1*par.Dsch12{p})*PP2';
            sumPP1 = PP1*ones(rr,1);
            tmp = sumPP1*sumPP1' + tmp2 + tmp2';
        else
            tmp = sparse(n,n);
        end
        vv = vectriu(pblk,tmp);
    elseif strcmp(pblk{1},'l')
        vv = par.Dsch12{p};
    end
    len = length(vv);
    D  = spdiags(vv,0,len,len);
    dd = dd + sum(At{p}.*(D*At{p}))';
end
%%*****************************************************************
%%*****************************************************************
%% vectriu: stack the upper triangular part of a matrix
%%          into a vector.
%%
%% xvec = vectriu(blk,x);
%%
%% SDPNAL:
%% Copyright (c) 2008 by
%% Xinyuan Zhao, Defeng Sun, and Kim-Chuan Toh
%%*****************************************************************

function xvec = vectriu(blk,x);

if ~iscell(x)
    r = blk{1,2};
    ismt = isempty(x);
    if strcmp(blk{1,1},'s') & (r > 0);
        r2 = r*(r+1)/2;
        xvec = zeros(r2,1);
        idx = 0;
        for j = 1:r
            xvec(idx+[1:j]) = x(1:j,j);
            idx = idx + j;
        end
    elseif strcmp(blk{p,1},'l') & (r > 0);
        xvec = x;
    end
else
    xvec = [];
    for p = 1:size(blk,1);
        xvec = [xvec; vectriu(blk(p,:),x{p})];
    end
end
%%*****************************************************************
