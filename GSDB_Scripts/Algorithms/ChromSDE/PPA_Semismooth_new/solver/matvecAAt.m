%%
%%  By = msk.*y + sig*AAt(y)
%%

function By = matvecAAt(At,par,y)

m = length(y);
By = zeros(m,1);
if (norm(y)==0);  return; end
%%
Aty = Atyfun(At,y,par.p,par.q);
By = par.msk.*y + par.sig*AXfun(At,Aty);

