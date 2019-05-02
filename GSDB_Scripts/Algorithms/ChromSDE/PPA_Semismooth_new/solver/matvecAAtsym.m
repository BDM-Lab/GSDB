%%************************************************************************
%% alternatingmatvec: matrix-vector multiply.
%% matrix B = tau*I + sig*A*At
%%************************************************************************

function By = matvecAAtsym(blk,At,par,y)

N = length(y);
By = zeros(N,1);
if (norm(y) == 0); return; end
%%
Aty = Atyfunsym(blk,At,y);
By  = par.sig*AXfunsym(blk,At,Aty);
%%
if isfield(par,'msk');
    By = By + par.msk.*y;
end
%%************************************************************************
