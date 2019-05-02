%%*****************************************************************
%% precondEDM: compute an approximation to
%%   (1/sigma)V approx = [M, q; q', alpha],
%%   q = A(H) and alpha = sum(sum(H));
%%   PrcondData provide the indices for generating the matrix J
%%   it is easier to generate Jt (the transpose of J)
%%   by default, we use the eigenvalue decomposition to
%%   compute the inverse of "Lambda + Jhat'*inv(D)*Jhat"
%%*****************************************************************

function L = precondEDM(blk,At,PrecondData,par)

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
        %%----------------------------------------------
        hbar = (PrecondData.ww).*h(PrecondData.indx_hbar);
        d = 1/sig + 2*hbar;
        dinv = 1./d;
        lenq = length(q);
        D = spdiags(dinv,0,lenq,lenq);
        %%
        hh = h(PrecondData.indx_h); 
        VV = (PrecondData.sqrtw).*sqrt(hh);
        row_indx = PrecondData.row_indx;
        col_indx = PrecondData.col_indx;
        Jt = spconvert([row_indx,col_indx,VV; n,lenq,0]);
        J = Jt';
        Dq = D*q; qDq = q'*Dq;
        if (alpha < 1e-1)
            precond_option = 2; fprintf('q1');
        elseif (abs(qDq-alpha) < 1e-1)
            precond_option = 2; fprintf('q2');            
        else
            precond_option = 1;
        end             
        if (precond_option ==1)
            Jhat = [J, q];
            DJhat = D*Jhat;
            Lambda = speye(n+1,n+1);
            Lambda(end,end) = -alpha;
            M = Jhat'*DJhat + Lambda;
        elseif (precond_option ==2)
            DJhat = D*J;
            M = Jt * DJhat + speye(n,n);
            fprintf('(II)');
        end
        inverse_option = 1;
        if (inverse_option == 1) 
            if issparse(M); M = full(M); end;
            L.matfct_options = 'lu';
            [L.L,L.U,L.perm] = lu(M,'vector');
            if (min(abs(diag(L.U))) < 1e-8)
               fprintf('switch to eigen-decomposition')
               inverse_option = 2;
            end
        end
        if (inverse_option == 2)
            % eigenvalue decomposition of M
            M = 0.5*(M+M');
            [P,eigM] = mexeig(full(M));
            dm = diag(eigM);
            fprintf(' | %2.1e, %2.1e, ',min(abs(full(dm))),max(abs(full(dm))));
            fprintf('%d',length(find(dm<0)));
            dm = 1./dm;
            L.P = P; L.Pt = P';
            L.dm = dm;
        end
        L.alpha = alpha;
        L.dinv  = dinv;
        L.DJhat = DJhat;
        L.DJhat_transpose = DJhat';
        L.precond_option = precond_option;
        L.inverse_option = inverse_option;
        if (precond_option == 1)
            invSq = invmatvecfun(L,q); 
            vec = zeros(lenq+1,1);
            vec(1:end-1) = -invSq/alpha;
            vec(end) = (1/alpha) + q'*invSq/alpha^2;
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
