%%********************************************************************
%%  PPASolver:
%%  Solve: min{ 0.5*norm(AX-b)^2 + rho*sum(svd(X)) : BX = d;  X\in R^{p*q}}
%%  Assume: p<= q
%%  AhatT = [AT  BT]  bhat = (b;d);  yhat =(y; z);
%%  [X,y,Z,info,runhist] =PPASolver(At,b,numeq,p,q,rho_target,OPTIONS,X,y,Z)
%%  numeq:  the number of extra constraints in BX = d.
%%*******************************************************************

function [obj,X,y,Z,runhist,W] =PPASolver(At,C,b,numeq,p,q,rho_target,OPTIONS,X,y,Z)

if (nargin < 8); OPTIONS = []; end
matlabversion = sscanf(version,'%f');
matlabversion = matlabversion(1);
warning off
randnstate = randn('state');
randn('state',0);

if (p>q)
    error(' By default: p should be less than or equal to q');
end

tol = 1e-6;
sig = 10;
maxiter       = 100;
maxitersub    = 10;
precond       = 1;
maxitpsqmr    = 200;
stagnate_check_psqmr = 0;
scale_data    = 2;
continuation  = 1;
printlevel    = 1;
use_QR = 0;  % use economical SVD
% plotyes       = 0;
if isfield(OPTIONS,'tol');        tol           = OPTIONS.tol; end
if isfield(OPTIONS,'sigma');      sig           = OPTIONS.sigma; end
if isfield(OPTIONS,'maxiter');    maxiter       = OPTIONS.maxiter; end
if isfield(OPTIONS,'maxitersub'); maxitersub    = OPTIONS.maxitersub; end
if isfield(OPTIONS,'precond');    precond       = OPTIONS.precond; end
if isfield(OPTIONS,'maxitpsqmr'); maxitpsqmr    = OPTIONS.maxitpsqmr; end
if isfield(OPTIONS,'stagnate_check_psqmr')
    stagnate_check_psqmr = OPTIONS.stagnate_check_psqmr;
end
if isfield(OPTIONS,'scale_data'); scale_data    = OPTIONS.scale_data; end
if isfield(OPTIONS,'continuation'); continuation= OPTIONS.continuation; end
if isfield(OPTIONS,'printlevel'); printlevel    = OPTIONS.printlevel; end
% if isfield(OPTIONS,'plotyes');    plotyes       = OPTIONS.plotyes; end
if isfield(OPTIONS,'use_QR');  use_QR = OPTIONS.use_QR; end

%%
tstart = clock;
if (printlevel)
    fprintf('\n size of matrix:  p= %2.0d, q=%2.0d', p,q);
    fprintf('\n num. of constraints = %2.0d',length(b));
end
m   = length(b);
b0  = b;   C0 = C;
msk = [ones(m-numeq,1); zeros(numeq,1)];
par.p =p;   par.q = q;
par.msk = msk;
par.precond = precond;
%%
rho_target_0 = rho_target;
if (continuation)
    % rho = max(rho_target,min(2,5e-2*normAtb));
    rho = 100*rho_target;
else
    rho = rho_target;
end
par.rho = rho;

%% scaling
%%
if (scale_data)
    % normA0 = min(1e12, max(1,norm(At,'fro')));
    % normb0 = max(1,norm(b));
    % normA0 = 1;
    % normb0 = min(1e12, max(1,sqrt(norm(b))));
    % normA0 = min(1e12, max(1,0.1*sqrt(norm(At,'fro'))));
    normA0 = min(1e12, max(1,sqrt(norm(At,'fro'))));  % this scaling works better
    normb0 = 1;
    normC0 = min(1e12, max(1,norm(C,'fro')));
    At = At/normA0;
    b = b/normb0;
    C = C/(normA0*normb0);    % C = C/normC0;
    par.rho =  par.rho/(normA0*normb0);   %  scaling ''par.rho'', not ''rho''
    objscale = normb0^2;
else
    normA0 = 1; normb0 = 1; normC0 = 1;
    objscale = 1;
end
normb = max(1,norm(b)); normC = max(1,norm(C,'fro'));
normA = max(1,norm(At,'fro'));
if (nargin < 11)
    X = zeros(p,q);  Z = zeros(p,q);  y = zeros(m,1);
    % normX = 0;
    generate_initial_point = 1;
else
    X = X*(normA0/normb0);
    svdX = mexsvd(X);
    nuclearX = sum(svdX);
    % normX = norm(svdX);
    rankX = length(find(svdX>0));
    y = y/normb0;
    Z = Z/(normA0*normb0);
    iter = 0;
    generate_initial_point = 0;
end
%%  use an alternating direction method to compute initial X,y,Z
% sig0 = 10;
sig0 = 10^3;
if (p>=300)
    maxiterAlt = 30;
elseif (p>= 200)
    maxiterAlt = 50;
else
    maxiterAlt = 50;
end
if (generate_initial_point)
    diagAAt = sum(At.*At)';
    diagAAt = max(1e-4,full(diagAAt));
    if (printlevel)
        fprintf('\n-----------------------------------------------------------------------------------------------------------------');
        fprintf('---------------');
        fprintf('\n iter     pinf            dinf            ratio      psqmr    rankX     sig');
        fprintf('\n-----------------------------------------------------------------------------------------------------------------');
        fprintf('---------------');
    end
    normX = norm(X,'fro');
    AX = AXfun(At,X);
    AC = AXfun(At,C);
    sigAZ = sig0*AXfun(At,Z);
    Rptmp = b - AX;
    for iter = 1: maxiterAlt
        par.sig = sig0;
        Xold = X; normXold = normX;
        yold = y;
        sigAZold = sigAZ;
        rhs = Rptmp + sig0*AC - sigAZ;
        diagH = msk + sig0*diagAAt;
        %% diagH = 1+ sig0*diagAAt;
        L.invdiagM = 1./diagH;
        [y,resnrm,solve_ok] = ...
            psqmr('matvecAAt',At,rhs,par,L,1e-3*norm(rhs),100,0);
        sigAty = Atyfun(At,sig0*y,p,q);
        W = X - sig0*C + sigAty;
        [X,svdX,rankX,sigZ] = threshold(W,par);
        AX = AXfun(At,X);
        normX = norm(svdX);
        nuclearX = sum(svdX);
        sigAZ = AXfun(At,sigZ);
        Rptmp = b - AX;
        Rp = Rptmp - msk.*y;
        normRp = norm(Rp);
        normRd = sqrt(normX^2 +normXold^2 - 2*sum(sum(X.*Xold)))/sig0;
        prim_infeas = normRp/normb;
        dual_infeas = normRd/normA;
        % ratio = prim_infeas/dual_infeas;
        ratio = prim_infeas/(eps+dual_infeas);
        if (printlevel)
            fprintf('\n %2.0d     %3.2e    %3.2e     %3.2e    %2.0d          %2.1d       %3.1e',...
                iter,prim_infeas,dual_infeas,ratio,length(resnrm),rankX,sig0);
        end
        if (rem(iter,3)==1)
            if (iter==1) & (ratio > 10)
                const = 10;
                X = Xold;  normX = normXold;
                y = yold;   sigAZ = sigAZold;
                Rptmp = b - AXfun(At,Xold);
                if (printlevel); fprintf('  ratio = %3.2e ',ratio); end
            else
                const = 2;
            end
            if (ratio < 0.1)
                sig0 = min(1e6,const*sig0);
                % sig0 = min(1e3,const*sig0);
            elseif (ratio > 10)
                sig0 = max(1e-2,sig0/const);
            end
        end
        if (max(prim_infeas,dual_infeas) < 1e-2)
            break;
        end
    end
end
%%
if (scale_data ==2)
    normX = norm(X,'fro');
    Z = sigZ/sig0;
    normZ = norm(Z,'fro');
    % normZ = norm(sigZ,'fro')/sig0;
    if (normZ > 1e-2)
        norm_ratio = full(normX/normZ);
        % norm_ratio = sqrt(norm_ratio);  % try this one later!
        % norm_ratio = 0.5*sqrt(norm_ratio);
        norm_ratio = sqrt(norm_ratio)/sqrt(2);
    else
        Aty = Atyfun(At,y,p,q);
        normtmp = max(1e-2,norm(C-Aty,'fro'));
        norm_ratio = full(normX/normtmp);
        % norm_ratio = sqrt(norm_ratio);
        % norm_ratio = 0.5*sqrt(norm_ratio);
        norm_ratio = sqrt(norm_ratio)/sqrt(2);
    end
    X = X/norm_ratio;  nuclearX = nuclearX/norm_ratio;
    Z = norm_ratio*Z;
    At = At*norm_ratio; C = C*norm_ratio;
    normA0 = normA0/norm_ratio;
    par.rho = par.rho * norm_ratio;
    normA = max(1,norm(At,'fro'));
end
%%
if (scale_data==2)
    sig = 1; % 5
    % sig = 10;
else
    norm_ratio = 1;
    normX = norm(X, 'fro');  normZ = norm(Z,'fro');
    sig = max(10,max(1e-2,normX)/max(1e-2,normZ));
end

% sig = sig/norm_ratio;  % check this one later!

AX   = AXfun(At,X);
Aty  = Atyfun(At,y,p,q);
Rp   = b-msk.*y-AX;
Rd   = C - Aty - Z;  % this equals to Rd in line 145
objtmp = 0.5*(norm(msk.*y))^2;
trCX = sum(sum(C.*X));
obj  = objscale*[objtmp+par.rho*nuclearX+trCX, -objtmp+b'*y];
% relgap   = abs(diff(obj))/max(1,mean(abs(obj)));
relgap   = (obj(1)-obj(2))/max(1,mean(abs(obj)));
prim_infeas  = norm(Rp)/normb;
dual_infeas  = norm(Rd, 'fro')/normA;
% ttime(1) = etime(clock,tstart);
runhist.prim_obj(1) = obj(1);
runhist.dual_obj(1) = obj(2);
runhist.gap(1)     = abs(diff(obj));
runhist.relgap(1) = relgap;
runhist.prim_infeas(1) = prim_infeas;
runhist.dual_infeas(1) = dual_infeas;
runhist.cputime(1) = etime(clock,tstart);
runhist.sigma(1)   = sig;
runhist.itersub(1) = iter;
runhist.rankX(1) = rankX;
runhist.findstep(1) = 0;
runhist.psqrm(1) = 0;

%%
fprintf('\n *****************************************************')
fprintf('********************************')
fprintf('\n  PPASolver For Nonsymmetric Matrix Problems')
fprintf('\n  rho_target = %3.2e',rho_target_0);
fprintf('\n  continuation = %1.0f',continuation)
fprintf('\n  scale_data = %1.1d, precond = %1.1d',scale_data,precond);
fprintf('\n  normC = %2.1e, normA = %2.1e, normb = %2.1e',...
    norm(C,'fro'), norm(At,'fro'), norm(b));
fprintf('\n  number of initial iterations = %2.0d',iter);
fprintf('\n  sig0 = %3.2e,  iter = %2.1d',sig0,iter);
fprintf('\n  norm_XZratio = %2.1e',norm_ratio);
fprintf('\n  orig-normX = %2.1e, orig-normZ = %2.1e',normX,normZ);
fprintf('\n  normX = %2.1e, normZ = %2.1e', normX/norm_ratio, normZ*norm_ratio);
fprintf('\n ****************************************************');
fprintf('******************************************')
fprintf('\n  it   prim_inf    dual_inf         prim_obj          dual_obj')
fprintf('           relgap      time     sigma      rho      rankX')
fprintf('\n *****************************************************')
fprintf('*****************************************')
fprintf('\n %3.1d  %3.1e    %3.1e  %- 8.7e %- 8.7e   %- 3.2e   %5.1f',...
    0,prim_infeas,dual_infeas,obj(1),obj(2),relgap,runhist.cputime(1));
fprintf('    %3.1e  %3.2e    %2.0d',sig,rho,rankX);
if (rho <= rho_target)&(max(prim_infeas,dual_infeas) < tol)
    if (printlevel)
        fprintf('\n--------------------------------------------------------')
        fprintf('\n max(prim_infeas,dual_infeas) < %3.2e',tol);
        fprintf('\n--------------------------------------------------------')
    end
    return;
end
%%
%% main
%%
L           = [];
msg         = '';
printsub    = printlevel;
% normX = normX/norm_ratio;
AAt = At.*At;
% diagAAt = sum(At.*At)';
diagAAt = sum(AAt)';
diagAAt = full(diagAAt);
%    if (printlevel)
%       fprintf('\n  mindiagAAt, maxdiagAAt = %3.2e, %3.2e',...
%                   min(diagAAt),max(diagAAt));
%    end
% diagAAt = max(1e-4,full(diagAAt));
diagAAt = max(1e-4,diagAAt);
%%
par.b   = b;
par.precond = precond;
par.use_QR = use_QR;
matvecfname = 'matvec';
par.project = 'threshold2';

%%
for iter = 1:maxiter
    % tstart = clock;
    par.rho = rho/(normA0*normb0);
    par.sig = sig;
    if (dual_infeas < 1e-5)
        maxitersub = max(maxitersub,50);
    elseif (dual_infeas < 1e-3)
        maxitersub = max(maxitersub,40);
    elseif (dual_infeas < 1e-1)
        maxitersub = max(maxitersub,30);
    end
    %%-------------------------------------------------------------------
    %%  Semismooth Newton method for inner subproblem
    %%-------------------------------------------------------------------
    if (prim_infeas < 0.2*dual_infeas)
        const = 0.5;
    elseif (prim_infeas > dual_infeas)
        const = 0.1;
    else
        const = 0.2;
    end
    tolsub = min(1,const*dual_infeas);
    %%--------------------------------
    %% compute new X, W
    %%--------------------------------
    Xold = X;
    % normXold = normX;
    normXold = norm(Xold,'fro');
    yold = y;
    XmsigC = X-sig*C;
    sigAty = Atyfun(At,sig*y,p,q);
    W   = XmsigC + sigAty;
    [X,par] = feval(par.project,W,par);
    normX = norm(par.dtp);  % norm(par.svdX)
    % rankX = length(par.posidx);
    %%
    clear subhist
    ysub = y;
    subhist.findstep(1) = 0;
    subhist.psqmr(1) = 0;
    % subhist.y(:,1) = ysub;
    for itersub = 1:maxitersub
        Ly     = -0.5*norm(msk.*ysub)^2 + b'*ysub -normX^2/(2*sig);
        Rpsub  = b - msk.*ysub-AXfun(At,X);
        GradLy = Rpsub;
        normRpsub = norm(Rpsub);
        normRdsub = sqrt(normX^2 +normXold^2 - 2*sum(sum(X.*Xold)))/sig;
        priminf_sub = normRpsub/normb;
        dualinf_sub = normRdsub/normA;
        const = 0.2;   % const = 0.5;
        tolsub = max(min(1,const*dualinf_sub), 1e-2*tol);
        subhist.priminf(itersub)  = priminf_sub;
        subhist.dualinf(itersub)  = dualinf_sub;
        subhist.Ly(itersub)       = Ly;
        if (printsub)
            fprintf('\n      %2.0d  %- 11.10e %3.2e %3.2e %3.1f',...
                itersub,Ly,priminf_sub,dualinf_sub,const);
        end
        %----------------------------------------------------------
        Ly_ratio = 1;  tolLy = 1e-8;
        if (itersub > 1)
            Ly_ratio = (Ly-subhist.Ly(itersub-1))/(abs(Ly)+eps);
        end
        %%---------------------------------------------
        %% check termination for inner problem
        %%----------------------------------------------
        if (priminf_sub < tolsub)&(itersub > 1)
            msg = 'good termination:';
            fprintf('\n       %s  ',msg);
            fprintf(' gradLy=%3.2e, priminf_sub=%3.2e, tolsub=%3.2e',...
                norm(GradLy), priminf_sub, tolsub);
            break;
        elseif (max(priminf_sub,dualinf_sub) < tol)
            msg = sprintf('max(prim_infeas,dual_infeas) < %3.2e',tol);
            fprintf('\n       %s',msg);
            break;
        elseif (abs(Ly_ratio) < tolLy) &(itersub>10)
            msg = sprintf('Ly_ratio < %3.1e',tolLy);
            fprintf('\n       %s',msg);
            break;
        end
        %-----------------------------------------------------------------
        %% compute Newton direction
        %% precond = 1, precondition by diag(A*diag(PI')*At)
        %%--------------------------------------------------------------
        diag_precond_options = 2;
        if (use_QR)
            diag_precond_options = 2;
        end
        if (diag_precond_options == 1)
            L = diagapprox(AAt,par);
        elseif (diag_precond_options ==2)
            diagH = msk + sig*diagAAt;
            % diagH = 1+sig*diagAAt;
            % diagH = max(1e-4, full(diagH));
            L.invdiagM = 1./diagH;
        end
        if (dual_infeas > 1e-3) | (itersub <= 5)
            maxitpsqmr = max(maxitpsqmr,200);
        elseif (dual_infeas > 1e-4)
            maxitpsqmr = max(maxitpsqmr,300);
        elseif (dual_infeas > 1e-5)
            maxitpsqmr = max(maxitpsqmr,400);
        elseif (dual_infeas > 5e-6)
            maxitpsqmr = max(maxitpsqmr,500);
        end
        par.minitpsqmr = 2;  % 2;
        if (dual_infeas > 1e-4)
            stagnate_check_psqmr = max(stagnate_check_psqmr,20);
        else
            stagnate_check_psqmr = max(stagnate_check_psqmr,50);
        end
        if (itersub > 3 & all(subhist.solve_ok(itersub-[3:-1]) <= -1)) ...
                & (dual_infeas < 5e-5)
            stagnate_check_psqmr = max(stagnate_check_psqmr,100);
        end
        par.stagnate_check_psqmr = stagnate_check_psqmr;
        if (itersub > 1)
            prim_ratio = priminf_sub/subhist.priminf(itersub-1);
            dual_ratio = dualinf_sub/subhist.dualinf(itersub-1);
        else
            prim_ratio = 0; dual_ratio = 0;
        end
        rhs = GradLy;
        % tolpsqmr = min(5e-2,0.1*norm(rhs));
        %% test the following part later!
        if (iter < 2) & (itersub < 5)
            tolpsqmr = min(1,0.1*norm(rhs));
        else
            tolpsqmr = min(5e-2,0.1*norm(rhs));
            % tolpsqmr = min(0.1,0.1*norm(rhs));
            % tolpsqmr = min(5e-3,0.1*norm(rhs));
        end
        %%
        const2 = 1.0;
        if (itersub>1) & (prim_ratio>0.5 | priminf_sub>0.1*subhist.priminf(1))
            const2 = 0.5*const2;
        end
        if (dual_ratio > 1.1); const2 = 0.5*const2; end
        tolpsqmr   = const2*tolpsqmr;
        [dy,resnrm,solve_ok] = ...
            psqmr(matvecfname,At,rhs,par,L,tolpsqmr,maxitpsqmr,printlevel);
        psqmriter = length(resnrm);
        if (printsub)
            fprintf('| %3.1e %3.1e %3.0d ',tolpsqmr,resnrm(end),psqmriter);
            fprintf(' %2.1f  %2.0d ',const2, length(par.dtp)); % length(par.svdX)
        end
        %% line search
        if (itersub <= 3) & (dual_infeas > 1e-4) | (iter <= 3)
            step_options = 1;
        else
            % step_options = 2;
            step_options = 1;
        end
        steptol = 1e-5;
        % steptol = 5e-4;
        [par,ysub,W,X,alpha,iterstep] = ...
            findstep(At,par,ysub,W,X,dy,steptol,step_options);
        normX = norm(par.dtp);
        % normX = norm(par.svdX);
        %
        subhist.solve_ok(itersub) = solve_ok;
        subhist.psqmr(itersub)    = psqmriter;
        subhist.findstep(itersub) = iterstep;
        % subhist.y(:,itersub+1)    = ysub;
        if (printsub)
            fprintf(' %3.2e %2.0f',alpha,iterstep);
            if (Ly_ratio < 0); fprintf('-'); end
        end
    end
    %%----- end of inner iteration -------------------
    
    %% test later!
    %       if (itersub == maxitersub)
    % 	     [priminfmin,idxmin] = min(subhist.priminf);
    % 	     idx  = find(subhist.priminf < 1.1*priminfmin);
    %    	     ysub = subhist.y(:,max(idx));  %% need store ysubs
    %          XmsigC = Xold- sig*C;
    %          sigAty = Atyfun(At,sig*ysub,p,q);
    %          W = XmsigC + sigAty;
    %          [X,par] = feval(par.project,W,par);
    %          % [X,svdX,rankX] = threshold(W,par);
    %          % normX = norm(svdX);
    %          % nuclearX = sum(svdX);
    %          normX = norm(par.dtp);
    %       end
    
    y   = ysub;
    rankX = length(par.posidx);
    nuclearX = sum(par.dtp);  % nuclear norm of X
    normX = norm(par.dtp);
    AX  = AXfun(At,X);
    Rp  = b-msk.*y-AX;
    normRp = norm(Rp);
    normRd = sqrt(normX^2 + normXold^2 - 2*sum(sum(X.*Xold)))/sig;
    prim_infeas = normRp/normb;
    dual_infeas = normRd/normA;
    objtmp  = 0.5*norm(msk.*y)^2;
    trCX = sum(sum(C.*X));
    obj = objscale*[objtmp+par.rho*nuclearX+ trCX, -objtmp+b'*y];
    % relgap  = abs(diff(obj))/max(1,mean(abs(obj)));
    relgap   = (obj(1)-obj(2))/max(1,mean(abs(obj)));
    %%
    runhist.prim_obj(iter+1) = obj(1);
    runhist.dual_obj(iter+1) = obj(2);
    runhist.gap(iter+1)     = abs(diff(obj));
    runhist.relgap(iter+1)   = relgap;
    runhist.prim_infeas(iter+1) = prim_infeas;
    runhist.dual_infeas(iter+1) = dual_infeas;
    runhist.cputime(iter+1) = etime(clock,tstart);
    runhist.sigma(iter+1)   = sig;
    runhist.itersub(iter+1) = itersub-1;
    runhist.findstep(iter+1) = sum(subhist.findstep);
    runhist.psqmr(iter+1)   = sum(subhist.psqmr);
    runhist.rankX(iter+1) = rankX;
    %%
    if (printlevel)
        fprintf('\n %3.0d  %3.1e    %3.1e  %- 8.7e %- 8.7e   %- 3.2e',...
            iter,prim_infeas,dual_infeas,obj(1),obj(2),relgap);
        fprintf('   %5.1f  ',runhist.cputime(iter+1));
        fprintf('   %3.1e  %3.2e    %2.1d',sig,rho,rankX);
    end
    %%
    %% check for termination
    %%
    idx = [max(iter-1,1):iter];
    prim_infeas_ratio = runhist.prim_infeas(idx+1)./runhist.prim_infeas(idx);
    dual_infeas_ratio = runhist.dual_infeas(iter+1)/runhist.dual_infeas(iter);
    if (rho <= rho_target)
        if (max(prim_infeas,dual_infeas) < tol)
            msg = sprintf('max(prim_infeas,dual_infeas) < %3.2e',tol);
            break;
        elseif ((dual_infeas < 0.5*tol) & all(prim_infeas_ratio > 0.7)) ...
                | ((prim_infeas_ratio(end) > 0.7) & all(runhist.dual_infeas(iter:iter+1) < 0.5*tol))
            msg = sprintf('lack of progress in prim_infeas');
            break;
        elseif (dual_infeas < 0.5*tol) & all(dual_infeas_ratio > 0.7)
            msg = sprintf('dual_infeas deteriorated');
            break;
        elseif (dual_infeas < 0.8*tol)
            msg = sprintf('dual_infeas < %3.2e',0.8*tol);
            break;
        end
    end
    %%
    %% reduce rho, modify sig
    %%
    if (continuation)
        if (max(prim_infeas,dual_infeas) < 10*tol)
            const = 0;
        elseif (max(prim_infeas,dual_infeas) < 100*tol)
            const = 0.2;
        else
            const = 0.3;
        end
        rho = max(rho_target,const*rho);
    end
    if (dual_infeas > 5*tol)
        sigtol = 0.7;  % try 0.5 later!
    else
        sigtol = 0.5;
    end
    if (dual_infeas_ratio > sigtol)
        const4 = 2;
        sig = min(1e8,const4*sig);
    end
end  %%  end of outer iteration

if (iter == maxiter); msg = 'maximum iteration reached'; end
Z = (X-W)/sig;
normXscaled = norm(X,'fro');
normyscaled = norm(y);
normZscaled = norm(Z,'fro');
if (scale_data)
    X = (normb0/normA0)*X;   y = normb0*y;
    Z = (normb0*normA0)*Z;
    Rp = b0-msk.*y-normA0*AXfun(At,X);
    Rd   = C0-Atyfun(At,normA0*y,p,q)-Z;
    prim_infeas = norm(Rp)/max(1,norm(b0));
    % At = normA0*At;
    dual_infeas = norm(Rd,'fro')/max(1,normA0*norm(At,'fro'));
end
rel_gap = (obj(1)-obj(2))/max(1,mean(abs(obj)));
Num_SVD_Decomp = runhist.itersub(1) + iter + sum(runhist.findstep(2:end));
fprintf('\n--------------------------------------------------------')
fprintf('------------------------')
fprintf('\n %s',msg);
fprintf('\n--------------------------------------------------------')
fprintf('------------------------')
fprintf('\n primal objval = %9.8e',obj(1));
fprintf('\n dual   objval = %9.8e',obj(2));
fprintf('\n relative gap  = %3.2e',rel_gap);
fprintf('\n prim_infeas   = %3.2e',prim_infeas);
fprintf('\n dual_infeas   = %3.2e',dual_infeas);
fprintf('\n CPU time      = %3.1f',runhist.cputime(end));
fprintf('\n number of subproblems = %3.1f',sum(runhist.itersub(2:end)));
fprintf('\n average number of psqmr step per subproblem = %3.1f',...
    sum(runhist.psqmr(2:end))/sum(runhist.itersub(2:end)));
fprintf('\n average number of linesearch per subproblem = %3.1f',...
    sum(runhist.findstep(2:end))/sum(runhist.itersub(2:end)));
fprintf('\n total number of singular value decomposition = %2.0d',Num_SVD_Decomp);
fprintf('\n scaled-norm(X) = %3.1e, scaled-norm(y) = %3.1e,',...
    normXscaled,normyscaled);
fprintf(' scaled-norm(Z) = %3.1e',normZscaled);
fprintf('\n norm(X) = %3.1e,        norm(y) = %3.1e,',...
    norm(X,'fro'),norm(y));
fprintf('        norm(Z) = %3.1e',norm(Z,'fro'));
fprintf('\n----------------------------------------------------------------');
fprintf('-----------------------------------------------------------');
fprintf('\n');
randn('state',randnstate);
%%********************************************************************
