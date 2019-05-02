%%*****************************************************************
%% precondsym: compute an approximation to
%%   (1/sigma)V approx = [M, q; q', alpha],
%%   q = A(H) and alpha = sum(sum(H));
%%   PrcondData provide the indices for generating the matrix J
%%   it is easier to generate Jt (the transpose of J)
%%   by default, we use the eigenvalue decomposition to
%%   compute the inverse of "Lambda + Jhat'*inv(D)*Jhat"
%%*****************************************************************

function L = precondsym(blk,At,PrecondData,par)

sig = par.sig;

for p = 1:size(blk,1)
    pblk = blk(p,:);
    n = sum(pblk{2});
    if strcmp(pblk{1},'s')
        PP1 = (par.P1{p}).*(par.P1{p});
        PP2 = (par.P2{p}).*(par.P2{p});
        rr = size(par.P1{p},2);
        if (rr == n)
            sumPP1 = PP1*ones(rr,1);
            tmp = sumPP1*sumPP1';
        elseif (rr > 0)
            tmp2 = (PP1*par.Dsch12{p})*PP2';
            sumPP1 = PP1*ones(rr,1);
            tmp = sumPP1*sumPP1' + tmp2 + tmp2';
        else
            tmp = sparse(n,n);
        end
        tmpX{1} = tmp;
        AX = AXfunsym(blk,At,tmpX);
        q = AX(1:end-1); 
        alpha = AX(end);
        h = vectriu(pblk,tmp);
        %%
        hbar = (PrecondData.ww).*h(PrecondData.indx_hbar);
        d = 1/sig + 2*hbar;
        dinv = 1./d;
        lenq = length(q);
        D = spdiags(dinv,0,lenq,lenq);
        %%
        hh = h(PrecondData.indx_h); hh = sqrt(hh);
        VV = (PrecondData.sqrtw).*hh;
        row_indx = PrecondData.row_indx;
        col_indx = PrecondData.col_indx;
        Jt = spconvert([row_indx,col_indx,VV; n,lenq,0]);
        J = Jt';
        %
        Dq = D*q; qDq = q'*Dq;
        if (norm(q)<1e-4)
            precond_option = 2;
            fprintf('q:');
        elseif (abs(qDq - alpha) < 1e-4)
            precond_option = 2;
            fprintf('q::');
        else
            precond_option = 1;
        end
        if (precond_option ==1)
            Jhat = [J, q];
            DJhat = D*Jhat;
            Lambda = speye(n+1,n+1);
            Lambda(end,end) = - alpha;
            M = Jhat'*DJhat + Lambda;
            lenM = n+1;
            inverse_option = 2;  % eigenvalue decomposition;           
            if (inverse_option == 1)
              %% check the smallest eigenvalue of M
              %% if min(eig(M))<0, then M is not positive definite
                eigsopts.disp  = 0;
                eigsopts.tol   = 1e-6; % 1e-3
                eigsopts.maxit = 200;
                M2 = 0.5*(M+M'); % just in case
                [~,eigM2,flag] = eigs(M2,3,'SA',eigsopts);
                if (eigM2(1) < 0) % M is not positive definite
                    % precond_option = 2;
                    fprintf('*');
                    inverse_option = 2; % eigenvalue decomposition
                end
            end
        end
        if (precond_option ==2)
            DJhat = D*J;
            M = Jt * DJhat + speye(n,n);
            lenM = n;
            inverse_option = 2;
            fprintf('(II)');
        end
        if (inverse_option == 1)
            % Cholesky factorization of M
            if (nnz(M) < 0.2*lenM*lenM)
                use_spchol=1;
            else
                use_spchol=0;
            end
            if (use_spchol)
                [L.R,L.p,L.perm] = chol(M,'vector');
                L.Rt = L.R';
                L.matfct_options = 'spcholmatlab';
            else
                if issparse(M); M = full(M); end;
                L.matfct_options = 'chol';
                L.perm = [1: lenM];
                [L.R,L.p] = chol(M);
            end
            fprintf('chol');
        elseif (inverse_option == 2)
            % eigenvalue decomposition of M
            M = 0.5*(M+M');
            [P,eigM] = mexeig(full(M));
            dm = diag(eigM);
            fprintf(' | %2.1e, %2.1e, ',min(abs(full(dm))),max(abs(full(dm))));
            fprintf('%d',length(find(dm<0)));
            dm = 1./dm;
            L.P = P; L.Pt =P';
            L.dm = dm;
            % fprintf('eig');
        end
        L.alpha = alpha;
        L.dinv = dinv;
        L.DJhat = DJhat;
        L.DJhat_transpose = DJhat';
        L.precond_option = precond_option;
        L.inverse_option  = inverse_option;
        if (precond_option == 1)
            qbar = (1/alpha)*q;
            sqbar = invmatvecfun(L,qbar); % S^{-1}*qbar;
            vec = zeros(lenq+1,1);
            vec(1:end-1) = - sqbar;
            vec(end) = (1/alpha) + qbar'*sqbar;
            L.vec = vec;
            fprintf('(I)');
        end
    end
end

%%*****************************************************************
%%*****************************************************************
%% vectriu: stack the upper triangular part of a matrix
%%          into a vector.
%%
%% xvec = vectriu(blk,x);

function xvec = vectriu(blk,x)

if ~iscell(x)
    r = blk{1,2};
    ismt = isempty(x);
    if strcmp(blk{1,1},'s') & (r > 0);
        r2 = r*(r+1)/2;
        xvec = zeros(r2,1);
        idx = 0;
        for j = 1:r
            xvec(idx+[1:j]) = x(1:j,j);
            idx = idx + j;
        end
    elseif strcmp(blk{p,1},'l') & (r > 0);
        xvec = x;
    end
else
    xvec = [];
    for p = 1:size(blk,1);
        xvec = [xvec; vectriu(blk(p,:),x{p})];
    end
end
%%*****************************************************************
