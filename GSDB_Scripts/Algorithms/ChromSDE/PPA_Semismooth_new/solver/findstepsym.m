%%************************************************************
%% findstep: strong Wolfe line search
%%
%%************************************************************

function  [par,y,W,Wp,alp,iter] = ...
                       findstepsym(blk,At,par,y0,W0,Wp0,dy,tol,options)

if ~exist('tol'); tol = 1e-4; end
if ~exist('options'); options = 1; end
printlevel = 0;
if isfield(par,'printlevel'); printlevel = par.printlevel; end

maxit = ceil(log(1/(tol+eps))/log(2));
c1 = 1e-4; c2 = 0.9;

b   = par.b;
sig = par.sig;
msk = par.msk;
%%
g0  = dy'*(b-msk.*y0-AXfunsym(blk,At,Wp0));
objtmp = 0.5*norm(msk.*y0)^2;
Ly0 = -objtmp + (b'*y0) - ops(Wp0,'norm')^2/(2*sig);
if (g0 <= 0)
    alp = 0; iter = 0;
    y = y0; W = W0; Wp = Wp0;
    if (printlevel)
        fprintf('\n Need an ascent direction, %2.1e  ',g0);
        return;
    end
end
%%
sigAtdy = Atyfunsym(blk,At,sig*dy);
alp = 1; alpconst = 0.5;
for iter =1:maxit
    if (iter==1);
        alp = 1; LB = 0; UB = 1;
    else
        alp = alpconst*(LB+UB);
    end
    y = y0 + alp*dy;
    W = ops(W0,'+',ops(alp,'*',sigAtdy));
    [Wp,par] = feval(par.project,blk,W,par);
    galp = dy'*(b-msk.*y-AXfunsym(blk,At,Wp));
    objtmp = 0.5*norm(msk.*y)^2;
    Ly   = -objtmp + (b'*y) - ops(Wp,'norm')^2/(2*sig);
    if (iter==1)
        gLB = g0; gUB = galp;
        if (sign(gLB)*sign(gUB) > 0)
            if (printlevel); fprintf('|'); end
            return;
        end
    end
    if (abs(galp) < c2*abs(g0)) & (Ly - Ly0 - c1*alp*g0 > -1e-8/max(1,abs(Ly0)))
        if (options==1) | ((options == 2) & (abs(galp) < tol))
            if (printlevel); fprintf(':'); end
            return;
        end
    end
    if (sign(galp)*sign(gUB) < 0)
        LB = alp; gLB = galp;
    elseif (sign(galp)*sign(gLB) < 0)
        UB = alp; gUB = galp;
    end
end
if (printlevel); fprintf('m'); end
%%************************************************************
