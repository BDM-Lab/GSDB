%%*************************************************************************
%% linsysolvefun: Solve H*x = b
%%
%% x = linsysolvefun(L,b)
%% where L contains the triangular factors of H. 
%%*************************************************************************
 
  function x = linsysolvefun(L,b)
 
  x = zeros(size(b)); 
  for k=1:size(b,2)
     if strcmp(L.matfct_options,'chol')
        x(L.perm,k) = L.R \ (b(L.perm,k)' / L.R)'; 
     elseif strcmp(L.matfct_options,'spchol')
        x(:,k) = bwblkslvfun(L, fwblkslvfun(L,b(:,k)) ./ L.d);
     elseif strcmp(L.matfct_options,'spcholmatlab')
        x(L.perm,k) = mexbwsolve(L.Rt,mexfwsolve(L.R,b(L.perm,k))); 
     elseif strcmp(L.matfct_options,'lu')
        x(:,k) = L.u \ (L.l \ b(L.perm,k));
     elseif strcmp(L.matfct_options,'splu')     
        x(L.perm,k) = L.q*( L.u \ (L.l \ (L.p*b(L.perm,k))));
     end
  end
%%*************************************************************************
