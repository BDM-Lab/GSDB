%%
%% Matrix CG For alternating Newton Step
%%

function [p,flag,relres,iterk] = Newton_matrixcg(blk, At, b,par,L, tol,maxit)
% Initializations
n=blk{1,2};
r = b;  %We take the initial guess x0=0 to save time in calculating A(x0) 
n2b =norm(b,'fro');    % norm of b
tolb = tol * n2b;  % relative tolerance 
p = zeros(n,n);
flag=1;
iterk =0;
relres=1000; %%% To give a big value on relres

  if isfield(par,'P');  par.Pt = ops(par.P,'transpose'); end
  if isfield(par,'P1'); par.P1t = ops(par.P1,'transpose'); end
  if isfield(par,'P2'); par.P2t = ops(par.P2,'transpose'); end
  if isfield(par,'Dsch12'); par.Dsch21 = ops(par.Dsch12,'transpose'); end

% Precondition 
z =r;
rz1 = sum(sum(r.*z)); %r'*z; 
rz2 = 1; 
d = z;
% CG iteration
for k = 1:maxit
   if k > 1
       beta = rz1/rz2;
       d = z + beta*d;
   end
   %%
       rhs = - matvecfun(blk,At,par,d);
       tolpsqmr =  min(5e-2,0.1*norm(rhs));
       maxitpsqmr=20;
       [dy,resnrm1,solve_ok] = ...
            psqmr(par.matvecfname,blk,At,rhs,par,L,tolpsqmr,maxitpsqmr);
       w = matfun(blk,At,par,dy,d);   % w = A(d)
   %-----------------------------------------------------------------------
   denom = sum(sum(d.*w)); %d'*w;
   iterk =k;
   relres = norm(r,'fro')/n2b;   %relative residue = norm(r) / norm(b)
   if denom <= 0 
       p = d/norm(d,'fro'); % d is not a descent direction
       p = (p+p')/2;
       break % exit
   else
       alpha = rz1/denom;
       p = p + alpha*d;
       r = r - alpha*w;
   end
   z = r;
   if norm(r,'fro') <= tolb % Exit if Hp=b solved within the relative tolerance
       iterk =k;
       relres = norm(r,'fro')/n2b;  %relative residue =norm(r) / norm(b)
       flag =0;
       p = (p+p')/2; 
       break
   end
   rz2 = rz1;
   rz1 = sum(sum(r.*z));  %r'*z;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%************************************************************************
%% For alternating Newton Step
%% matrix B = DPi*H
%%************************************************************************

   function B = matfun(blk,At,par,dy,H)
   
   N=blk{1,2};  %length(H);
   if (norm(H,'fro') == 0) & (norm(dy)==0)
       B = zeros(N,N);
       return;
   end
 %%
   B = zeros(N,N);
   for p = 1:size(blk,1)
      pblk = blk(p,:);
      n = sum(pblk{2});
      if strcmp(pblk{1},'s')
         rr = size(par.P1{p},2);
         sigAtdy = Atyfun(pblk,At{p},par.sig*dy);
         Aty = H + sigAtdy;
         if (rr > 0 & rr < n) 
            if (rr <= n/2) 
               tmp0 = par.P1t{p}*Aty;
               tmp1 = (tmp0*par.P1{p})*par.P1t{p};         
               tmp2 = par.Dsch12{p}.*(tmp0*par.P2{p});
               tmp2 = tmp2*par.P2t{p}; 
               tmp3 = par.P1{p}*(0.5*tmp1 + tmp2);
               B = B + (tmp3+tmp3');
            else
               tmp0 = par.P2t{p}*Aty;
               tmp1 = (tmp0*par.P2{p})*par.P2t{p};         
               tmp2 = (1-par.Dsch21{p}).*(tmp0*par.P1{p});
               tmp2 = tmp2*par.P1t{p}; 
               tmp3 = par.P2{p}*(0.5*tmp1 + tmp2);
               B = B + (Aty-tmp3-tmp3');
            end
         elseif (rr == n)
               B = B + Aty;
         end
      elseif strcmp(pblk{1},'l')   % need modification
	     if (~isempty(par.Dsch12{p}))
            tmp = par.Dsch12{p}.*(At{p}*y); 
 	        By = By + par.sig*(tmp'*At{p})';  
         end
      end
   end
   
   B = (H-B)/par.sig;

%%************************************************************************

%%************************************************************************
%% For alternating Newton Step
%% matvecfun: matrix-vector multiply.  
%% matrix B = A*DPi*At
%%                 =A(PxP) Dsch (PtxPt)At 
%%************************************************************************

   function By = matvecfun(blk,At,par,H)
   
   N = par.m;
   if (norm(H,'fro') == 0); By = zeros(N,1); return; end
 %%
   By = zeros(N,1);
   for p = 1:size(blk,1)
      pblk = blk(p,:);
      n = sum(pblk{2}); 
      if strcmp(pblk{1},'s')
         rr = size(par.P1{p},2);
         Aty = H;   % Atyfun(pblk,At{p},y);
         if (rr > 0 & rr < n) 
            if (rr <= n/2) 
               tmp0 = par.P1t{p}*Aty;
               tmp1 = (tmp0*par.P1{p})*par.P1t{p};         
               tmp2 = par.Dsch12{p}.*(tmp0*par.P2{p});
               tmp2 = tmp2*par.P2t{p}; 
               tmp3 = par.P1{p}*(0.5*tmp1 + tmp2);
               By = By + AXfun(pblk,At{p},tmp3+tmp3');
            else
               tmp0 = par.P2t{p}*Aty;
               tmp1 = (tmp0*par.P2{p})*par.P2t{p};         
               tmp2 = (1-par.Dsch21{p}).*(tmp0*par.P1{p});
               tmp2 = tmp2*par.P1t{p}; 
               tmp3 = par.P2{p}*(0.5*tmp1 + tmp2);
               By = By + AXfun(pblk,At{p},Aty-tmp3-tmp3');
            end
         elseif (rr == n)
               By = By + AXfun(pblk,At{p},Aty);
         end
      elseif strcmp(pblk{1},'l')       % need modification
	     if (~isempty(par.Dsch12{p}))
            tmp = par.Dsch12{p}.*(At{p}*y); 
 	        By = By + par.sig*(tmp'*At{p})'; 
         end
      end
   end

%%************************************************************************