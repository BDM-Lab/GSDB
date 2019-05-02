%%************************************************************************
%% For SVD use nonsymmetric case
%% matvec: matrix-vector multiply.
%% matrix B = diag(msk) + sig* (A*DP_rhosigma*At)
%%************************************************************************

function By = matvec(At,par,y)

N = length(y);
p = par.p;   q = par.q;
if (norm(y) == 0); By = zeros(N,1); return; end
%%
By = zeros(N,1);
r = length(par.posidx);
if (r==0); By = par.msk.*y; return; end
%%
Aty = Atyfun(At,y,p,q);
if (r>0)&(r<p)
    Ur = par.U(:,[1:r]);    Urt = Ur';
    Us= par.U(:,[r+1:p]);
    V1r = par.V1(:,[1:r]);     V1s = par.V1(:,[r+1:p]);
    V1rt = par.V1t([1:r],:);
    % V1st = par.V1t([r+1:p],:);
    %%
    tmp0 = Aty*V1r;   tmp1 = Aty*V1s;
    Hrr = Urt*tmp0;   Hrs = Urt*tmp1;
    Hsr = Us'*tmp0;
    % Hrrt = Hrr';
    M11 = Hrr + Hrr' + par.Gamma11.*(Hrr - Hrr');
    Omega12 = par.Omega12;    Omega12t = par.Omega12t;
    Gamma12 = par.Gamma12;  Gamma12t = par.Gamma12t;
    M12 = (Omega12 + Gamma12).*Hrs;
    M12 = M12 + (Omega12 - Gamma12).*Hsr';
    M21 = (Omega12t + Gamma12t).*Hsr;
    M21 = M21 + (Omega12t - Gamma12t).*Hrs';
    tmpM = [M11 M12];
    tmp2 = tmpM*par.V1t;
    tmp2 = Ur*tmp2;
    if (r<=p/2)
        tmp3 = Us*M21;
        tmp3 = tmp3*V1rt;
    else
        tmp3 = M21*V1rt;
        tmp3 = Us*tmp3;
    end
    DPAty = 0.5*(tmp2+tmp3);
    if isfield(par,'HouseQ')
        tmp4 = mexAV2(Aty,par.HouseQ);
        tmp4 = Urt*tmp4;
        tmp5 = par.Beta.*tmp4;
        tmp5 = mexAV2(tmp5, par.HouseQ, 1);
        DPAty = DPAty + Ur*tmp5;
    elseif isfield(par,'V2')
        tmp4 = Aty*par.V2;
        tmp4 = Urt*tmp4;
        tmp5 = par.Beta.*tmp4;
        tmp5 = tmp5*par.V2t;
        DPAty = DPAty + Ur*tmp5;
    end
    By = par.msk.*y + par.sig*AXfun(At,DPAty);
elseif (r == p)
    tmp0 = Aty*par.V1;
    % Ut = (par.U)';
    H = (par.U)'*tmp0;     Ht = H';
    Hsym = 0.5*(H + Ht);
    Hskew = 0.5*(H - Ht);
    tmp1 = par.Gamma.*Hskew;
    tmp1 = Hsym + tmp1;
    tmp2 = par.U*tmp1;
    DPAty = tmp2*par.V1t;
    if isfield(par,'HouseQ')
        tmp3 = mexAV2(Aty, par.HouseQ);
        tmp3 = (par.U)'*tmp3;
        tmp3 = par.Beta.*tmp3;
        if ((q-p)<p)
            tmp3 = par.U*tmp3;
            DPAty = DPAty + mexAV2(tmp3, par.HouseQ,1);
        else
            tmp3 = mexAV2(tmp3, par.HouseQ,1);
            DPAty = DPAty + par.U*tmp3;
        end
    elseif isfield(par,'V2')
        tmp3 = Aty*par.V2;
        tmp3 = (par.U)'*tmp3;
        tmp3 = par.Beta.*tmp3;
        if ((q-p)<p)
            tmp3 = par.U*tmp3;
            DPAty = DPAty + tmp3*par.V2t;
        else
            tmp3 = tmp3*par.V2t;
            DPAty = DPAty + par.U*tmp3;
        end
    end
    By = par.msk.*y + par.sig*AXfun(At,DPAty);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


