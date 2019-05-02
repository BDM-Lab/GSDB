%%
%% soft-thresholding operator:
%%

function [Wp,svdWp,rankWp,Wn] = threshold(W,par)

tmpversion = version('-release');
matlabrelease = str2num(tmpversion(1:4));

[p,q] = size(W);
postol = 1e-12;  % postol = 1e-16;
if isfield(par,'sig');
    rho = par.rho*par.sig;  % rho = rho*sigma;
else
    rho = par.rho;
end

if (matlabrelease < 2011)
    [U,S,V] = mexsvd(W);  % economical SVD by default
    s = diag(S);
    %%
    %%  The lapack function DGESDD could return negative or NaN singular values
    %%
    z = randn(q,1);
    tmpz = U*(s.*(V'*z));
    Wz = W*z;  normWz = norm(Wz);
    ratio = norm(tmpz - Wz)/(1+normWz);
    if isnan(ratio) || (ratio > 1e-12)
        [U,S,V] = svd(X,'econ');
        if (~isnan(ratio))
            fprintf(' mexsvd ratio = %3.2e ', ratio);
        else
            fprintf(' N ');
        end
    end
    s = diag(S);
else
    [U,S,V] = svd(W,'econ');
    s = diag(S);
end

% s = diag(S);
stmp = s - rho;
idx = find(stmp> postol);
len = length(idx);
if (len > 0)
    UU = U(:,idx);
    VV = V(:,idx);
    svdWp = stmp(idx);
    SS = spdiags(svdWp,0,len,len);
    Wp = UU*SS*VV';
    if (nargout==4)
        Wn = Wp - W;
    end
else
    Wp = zeros(p,q);
    svdWp = [];
    if (nargout==4)
        Wn = - W;
    end
end
rankWp = len;

