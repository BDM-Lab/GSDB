%%*************************************************************************
%% psqmr:  preconditioned symmetric QMR with left (symmetric) preconditioner.
%%
%% b = rhs vector.
%% resnrm = norm of qmr-generated residual vector b-Ax.
%%

%% M=M1M2.  Here M2 = I , left preconditioning
%%*************************************************************************

function  [x,resnrm,solve_ok] = psqmrsym(matvecfname,blk,At,b,par,L,tol,maxit,x0)

N = length(b);
if ~exist('maxit'); maxit = max(50,sqrt(length(b))); end;
if ~exist('tol'); tol = 1e-6*norm(b); end;
if ~exist('x0'); x0 = zeros(N,1); end;
if ~exist('L');  L = []; end;

solve_ok = 1;
stagnate_check = 20;
miniter = 2;
if isfield(par,'minitpsqmr'); miniter = par.minitpsqmr; end
if isfield(par,'stagnate_check_psqmr');
    stagnate_check = par.stagnate_check_psqmr;
end
%%
if isfield(par,'P');  par.Pt = ops(par.P,'transpose'); end
if isfield(par,'P1'); par.P1t = ops(par.P1,'transpose'); end
if isfield(par,'P2'); par.P2t = ops(par.P2,'transpose'); end
if isfield(par,'Dsch12'); par.Dsch21 = ops(par.Dsch12,'transpose'); end
%%
x = x0;
if (norm(x) > 0)
    Aq = feval(matvecfname,blk,At,par,x);
else
    Aq = zeros(N,1);
end
r = b-Aq;
err = norm(r); resnrm(1) = err; minres = err;
%%
if isnumeric(par.precond)
    q = precondfun(blk,At,par,L,r);
else
    q = feval(par.precond.fname,r,par.precond);
end
tau_old  = norm(q);
rho_old  = r'*q;
theta_old = 0;
d = zeros(N,1);
res = r; Ad = zeros(N,1);
%%
%% main loop
%%
tiny = 1e-30;
for iter = 1:maxit
    Aq = feval(matvecfname,blk,At,par,q);
    sigma = q'*Aq;
    if (abs(sigma) < tiny)
        solve_ok = 2;
        %%fprintf('sm');
        break;
    else
        alpha = rho_old/sigma;
        r = r - alpha*Aq;
    end
    if isnumeric(par.precond)
        u = precondfun(blk,At,par,L,r);  %% M2 =I, left preconditioning
    else
        u = feval(par.precond.fname,r,par.precond);
    end
    %%
    theta = norm(u)/tau_old; c = 1/sqrt(1+theta^2);
    tau = tau_old*theta*c;
    gam = (c^2*theta_old^2); eta = (c^2*alpha);
    d = gam*d + eta*q;
    x = x + d;
    %%----- stopping conditions ----
    Ad = gam*Ad + eta*Aq;
    res = res - Ad;
    err = norm(res); resnrm(iter+1) = err;
    if (err < minres); minres = err; end
    if (strcmp(matvecfname,'matvecAAtsym'))
        if (err < tol) & (iter > miniter); break; end
    elseif(strcmp(matvecfname,'matvecsym'))
        if (err < tol) & (iter > miniter) & (b'*x > 0); break; end
    end
    if (iter > stagnate_check) & (iter > 10)
        ratio = resnrm(iter-9:iter+1)./resnrm(iter-10:iter);
        if (min(ratio) > 0.997) & (max(ratio) < 1.003)
            solve_ok = -1; fprintf('s')
            break;
        end
    end
    %%-----------------------------
    if (abs(rho_old) < tiny)
        solve_ok = 2;
        %%fprintf('sm');
        break;
    else
        rho  = r'*u;
        beta = rho/rho_old;
        q = u + beta*q;
    end
    rho_old = rho;
    tau_old = tau;
    theta_old = theta;
end
if (iter == maxit); solve_ok = -2; end
%%
%%************************************************************************

function  q = precondfun(blk,At,par,L,r)

precond = 0;
if isfield(par,'precond'); precond = par.precond; end

if (precond == 0)
    q = r;
elseif (precond == 1)
    q = L.invdiagM.*r;
elseif (precond == 2)
    if strcmp(L.matfct_options,'chol')
        q(L.perm,1) = mextriang(L.L, mextriang(L.L,r(L.perm),2) ,1);
    elseif strcmp(L.matfct_options,'spchol')
        q = bwblkslvfun(L, fwblkslvfun(L,r) ./ L.d);
    elseif strcmp(L.matfct_options,'spcholmatlab')
        q(L.perm,1) = mexbwsolve(L.Rt,mexfwsolve(L.R,r(L.perm,1)));
    end
elseif (precond == 3)
    %% for RKE and EDM problems
    x = r(1:end-1); y = r(end);
    if (L.precond_option == 1)
        %% invPrecondMat = [invS, -invSq/alp; -invSq'/alp, 1/alp+q'*invSq/alp^2]
        %% L.vec(1:end-1) = -invSq/alp;
        %% L.vec(end) = 1/alp + q'*invSq/alp^2
        sx = invmatvecfun(L,x); % inv(S)*x;
        tmp = [sx; x'*L.vec(1:end-1)];
        q = tmp + y*L.vec;
    elseif (L.precond_option == 2)
        % invPrecondMat = [invM; 1]
        mx = invmatvecfun(L,x); % mx = inv(M)*x;
        q = [mx; y]; 
    end
    q = (1/par.sig) *q;
end
%%************************************************************************
