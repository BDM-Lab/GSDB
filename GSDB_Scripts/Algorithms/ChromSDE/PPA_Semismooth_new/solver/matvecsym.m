%%************************************************************************
%% matvec: matrix-vector multiply.
%% matrix B = diag(msk) + sig* (A*DPi*At)
%%          = diag(msk) + sig * A(PxP) Dsch (PtxPt)At
%%************************************************************************

function By = matvecsym(blk,At,par,y)

N = length(y);
if (norm(y) == 0); By = zeros(N,1); return; end
%%
By = zeros(N,1);
for p = 1:size(blk,1)
    pblk = blk(p,:);
    n = sum(pblk{2});
    if strcmp(pblk{1},'s')
        rr = size(par.P1{p},2);
        Aty = Atyfunsym(pblk,At{p},y);
        if (rr > 0 & rr < n)
            if (rr <= n/2)
                tmp0 = par.P1t{p}*Aty;
                tmp1 = (tmp0*par.P1{p})*par.P1t{p};
                tmp2 = par.Dsch12{p}.*(tmp0*par.P2{p});
                tmp2 = tmp2*par.P2t{p};
                tmp3 = par.P1{p}*(0.5*tmp1 + tmp2);
                By = By + par.sig*AXfunsym(pblk,At{p},tmp3+tmp3');
            else
                tmp0 = par.P2t{p}*Aty;
                tmp1 = (tmp0*par.P2{p})*par.P2t{p};
                tmp2 = (1-par.Dsch21{p}).*(tmp0*par.P1{p});
                tmp2 = tmp2*par.P1t{p};
                tmp3 = par.P2{p}*(0.5*tmp1 + tmp2);
                By = By + par.sig*AXfunsym(pblk,At{p},Aty-tmp3-tmp3');
            end
        elseif (rr == n)
            By = By + par.sig*AXfunsym(pblk,At{p},Aty);
        end
    elseif strcmp(pblk{1},'l')
        if (~isempty(par.Dsch12{p}))
            tmp = par.Dsch12{p}.*(At{p}*y);
            By = By + par.sig*(tmp'*At{p})';
        end
    end
end
%%
if isfield(par,'msk')
    By = By + par.msk.*y;
end
%%************************************************************************
