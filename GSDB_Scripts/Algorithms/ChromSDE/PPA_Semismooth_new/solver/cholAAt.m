%%***********************************************************
%% cholAAt: compute Cholesky factorization of A*At. 
%%
%% L = Cholesky factor of AAt + sig*I
%%
%%***********************************************************

   function [L,AAt] = cholAAt(blk,At,m,sig); 

   if (nargin <= 3); sig = 0; end

   cachesize = 256; 
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1); 

   AAt = sparse(m,m);  
   for p = 1:size(blk,1)
      AAt = AAt + At{p,1}'*At{p,1}; 
   end
   pertdiag = 1e-13*ones(m,1); 
   AAt = AAt + spdiags(pertdiag,0,m,m);   
%%
   HH = AAt + sig*speye(m,m); 
%%
   if (nnz(HH) < 0.2*m*m); use_spchol=1; else; use_spchol=0; end
   if (use_spchol)
      if (matlabversion < 7.3)
         [Lsymb,flag] = symbcholfun(HH,cachesize);
         if (flag); existspcholsymb = 0; else; existspcholsymb = 1; end
         L = sparcholfun(Lsymb,HH);
         L.matfct_options = 'spchol';
      else
         [L.R,L.p,L.perm] = chol(HH,'vector'); 
         L.Rt = L.R'; 
         L.matfct_options = 'spcholmatlab'; 
      end
   else
      if issparse(HH); HH = full(HH); end;
      L.matfct_options = 'chol';
      L.perm = [1:m]; 
      [L.R,indef] = chol(HH);
   end
%%***********************************************************
