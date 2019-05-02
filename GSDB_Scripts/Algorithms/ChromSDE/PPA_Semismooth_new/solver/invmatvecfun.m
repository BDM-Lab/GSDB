%%******************************************************************
%% computer y = inv(S)*x;
%% inv(S) = invD - invDJhat*Ginv*invDJhat'
%%
%% For RKE and EDM problems
%%******************************************************************
function y = invmatvecfun(L,x)

if (norm(x) == 0); y = zeros(length(x),1); return; end

y = (L.dinv).*x;
r = L.DJhat_transpose*x;
if (L.inverse_option == 1)
    if strcmp(L.matfct_options,'chol')
       q(L.perm,1) = mextriang(L.R, mextriang(L.R,r(L.perm),2) ,1);
    elseif strcmp(L.matfct_options,'spcholmatlab')
       q(L.perm,1) = mexbwsolve(L.Rt,mexfwsolve(L.R,r(L.perm,1)));
    elseif strcmp(L.matfct_options,'lu')
       q = L.U\ (L.L\ r(L.perm)); 
    end
elseif (L.inverse_option == 2)
    q = L.Pt*r;
    q = (L.dm).*q;
    q = (L.P)*q;
end
y = y - L.DJhat*q;
%%******************************************************************