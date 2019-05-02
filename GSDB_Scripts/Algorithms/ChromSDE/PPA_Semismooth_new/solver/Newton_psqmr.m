%%*************************************************************************
%% Newton_psqmr:  preconditioned symmetric QMR with left (symmetric) preconditioner. 
%%
%% b = rhs Matrix.     x is a matrix
%% resnrm = norm of qmr-generated residual vector b-AX. 
%% subfuntions:   1.  Aq =matfun(blk,At,par,dy,q)
%%                          2.  A DPiH = matvecfun(blk,At,par,H)
%%
%%*************************************************************************

   function  [x,resnrm,solve_ok] = Newton_psqmr(blk,At,b,par,L,tol,maxit,x0) 

   N = blk{1,2}; 
   if ~exist('maxit'); maxit = 20; end;
   if ~exist('tol'); tol = 1e-2; end; 
   if ~exist('x0'); x0 = zeros(N,N); end;
   if ~exist('L');  L = []; end;

   tolb = tol* norm(b,'fro');
   solve_ok = 1; 
   stagnate_check = 20; 
   miniter = 2;
   if isfield(par,'minitpsqmr'); miniter = par.minitpsqmr; end
%%
   if isfield(par,'P');  par.Pt = ops(par.P,'transpose'); end
   if isfield(par,'P1'); par.P1t = ops(par.P1,'transpose'); end
   if isfield(par,'P2'); par.P2t = ops(par.P2,'transpose'); end
   if isfield(par,'Dsch12'); par.Dsch21 = ops(par.Dsch12,'transpose'); end
%%
   x = x0; 
%    if (norm(x) > 0) 
%       Aq = feval(matvecfname,blk,At,par,x);   
%    else
%       Aq = zeros(N,1);  
%    end
   r = b;  %b-Aq;
   err = norm(r,'fro'); resnrm(1) = err; minres = err; 
%%
%    if isnumeric(par.precond)
%       q = precondfun(blk,At,par,L,r); 
%    else
%       q = feval(par.precond.fname,r,par.precond); 
%    end
   q=r;
   tau_old  = norm(q,'fro');      
   rho_old  = sum(sum(r.*q));  % r'*q; 
   theta_old = 0; 
   d = zeros(N,N); 
   res = r; Ad = zeros(N,N);
%%      
%% main loop
%%
   tiny = 1e-30; 
   for iter = 1:maxit
       rhs = - matvecfun(blk,At,par,q);
       tolpsqmr =  min(5e-2,0.1*norm(rhs));
       maxitpsqmr=20;
       [dy,subresnrm,subsolve_ok] = ...
              psqmrsym(par.matvecfname,blk,At,rhs,par,L,tolpsqmr,maxitpsqmr);
      %Aq = feval(matvecfname,blk,At,par,dy,q);
      Aq = matfun(blk,At,par,dy,q);
      %---------------------------------------------------
       sigma = sum(sum(q.*Aq));  % q'*Aq; 
       if (abs(sigma) < tiny)
          solve_ok = 2; 
          %%fprintf('sm'); 
          x = (x+ x')/2;
          break;
       else
          alpha = rho_old/sigma; 
          r = r - alpha*Aq;
       end
%        if isnumeric(par.precond)
%           u = precondfun(blk,At,par,L,r); 
%        else
%           u = feval(par.precond.fname,r,par.precond); 
%        end
       u = r;
       %%
       theta = norm(u,'fro')/tau_old; c = 1/sqrt(1+theta^2); 
       tau = tau_old*theta*c;
       gam = (c^2*theta_old^2); eta = (c^2*alpha); 
       d = gam*d + eta*q;
       x = x + d; 
       %%----- stopping conditions ----
       Ad = gam*Ad + eta*Aq;
       res = res - Ad; 
       err = norm(res,'fro'); resnrm(iter+1) = err; 
       if (err < minres); minres = err; end
       if (err < tolb) & (iter > miniter) & (sum(sum(b.*x)) > 0)  % (b'*x > 0)
           x = (x+ x')/2;
           break;
       end  
       if (iter > stagnate_check) & (iter > 10)
          ratio = resnrm(iter-9:iter+1)./resnrm(iter-10:iter); 
          if (min(ratio) > 0.997) & (max(ratio) < 1.003)
             solve_ok = -1; fprintf('s')
             x = (x+ x')/2;
             break;
          end       
       end
       %%----------------------------- 
       if (abs(rho_old) < tiny)
          solve_ok = 2; 
          x = (x+ x')/2;
          %%fprintf('sm');
          break;
       else
          rho  = sum(sum(r.*u)); % r'*u; 
          beta = rho/rho_old; 
          q = u + beta*q; 
       end
       rho_old = rho; 
       tau_old = tau; 
       theta_old = theta; 
   end
   if (iter == maxit); solve_ok = -2; end
   
%  if (solve_ok ~= -1); fprintf(' '); end

%%************************************************************************

%%************************************************************************
%% For alternating Newton Step
%% matrix B = DPi *H
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
         sigAtdy = Atyfunsym(pblk,At{p},par.sig*dy);
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
      elseif strcmp(pblk{1},'l') % need modification
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
         Aty = H;
         if (rr > 0 & rr < n) 
            if (rr <= n/2) 
               tmp0 = par.P1t{p}*Aty;
               tmp1 = (tmp0*par.P1{p})*par.P1t{p};         
               tmp2 = par.Dsch12{p}.*(tmp0*par.P2{p});
               tmp2 = tmp2*par.P2t{p}; 
               tmp3 = par.P1{p}*(0.5*tmp1 + tmp2);
               By = By + AXfunsym(pblk,At{p},tmp3+tmp3');
            else
               tmp0 = par.P2t{p}*Aty;
               tmp1 = (tmp0*par.P2{p})*par.P2t{p};         
               tmp2 = (1-par.Dsch21{p}).*(tmp0*par.P1{p});
               tmp2 = tmp2*par.P1t{p}; 
               tmp3 = par.P2{p}*(0.5*tmp1 + tmp2);
               By = By + AXfunsym(pblk,At{p},Aty-tmp3-tmp3');
            end
         elseif (rr == n)
               By = By + AXfunsym(pblk,At{p},Aty);
         end
      elseif strcmp(pblk{1},'l')  % need modification
	     if (~isempty(par.Dsch12{p}))
            tmp = par.Dsch12{p}.*(At{p}*y); 
 	        By = By + par.sig*(tmp'*At{p})'; 
         end
      end
   end

%%************************************************************************


