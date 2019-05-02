%%******************************************************
%% refinedistances
%%
%% L, U must be upper-triangular distance matrices.
%%
%% Aest = refinedistances()
%%******************************************************

function [Aest,info] = refinepositions_new(Aorg,D,beta_param,maxIter,tolerance)

if ~exist('maxIter'); maxIter = 2000; end
if ~exist('tolerance'); tolerance = 1e-9; end
if ~exist('beta_param'); beta_param = 0.03; end;  %% best: 0.1-0.2

fid = 1;
D = triu(D);
idxall = [1:length(D)];
idxbad = [];
%%
for pass=[1:5]
    idxrun = setdiff(idxall,idxbad);
    Dtmp    = D(idxrun,idxrun);
    nAtoms  = size(Dtmp,1);
    nEdges  = nnz(Dtmp);
    [i,j,d] = find(Dtmp);
    eij     = sparse(i,1:nEdges,1,nAtoms,nEdges) ...
        - sparse(j,1:nEdges,1,nAtoms,nEdges);
    ee = ones(nAtoms,1)/sqrt(nAtoms);
    Astart = Aorg(:,idxrun)-(Aorg(:,idxrun)*ee)*ee';
    objstart = objfunc(Astart,d,eij,0);
    beta = beta_param*objstart/norm(Astart,'fro')^2;
    solveOK = 1;
    Aold = Astart;
    objOld = objfunc(Aold,d,eij,beta);
    Anew = Aold;
    for iter = 1:maxIter
        objNew = inf;
        Agrad  = objgrad(Aold,d,eij,beta);
        alpha  = 0.3;
        iBack  = 0;
        while (objNew > objOld) && (iBack < 20)
            %% Backtracking line search.
            iBack  = iBack + 1;
            alpha  = alpha / 2;
            Anew   = Aold - alpha * Agrad;
            objNew = objfunc(Anew,d,eij,beta);
        end
        if (abs(objOld-objNew)/(1+abs(objOld)) < tolerance)
            break
        end
        Aold   = Anew;
        objOld = objNew;
    end
    Anew = Anew-(Anew*ee)*ee';
    objend = objfunc(Anew,d,eij,0); %% 0 to exclude regularization term
    expansion_ratio = sqrt(sum(Anew.*Anew))./sqrt(sum(Astart.*Astart));
    expansion_idx = find(expansion_ratio > 3.0);
%     if (~isempty(expansion_idx))
%         Anew = Astart;
%         solveOK = 0;
%         fprintf(fid,'***** refined solution has expanded too much, discarded');
%         fprintf(fid,':%3.2e, %2.0f\n',max(expansion_ratio),length(expansion_idx));
%         if (length(expansion_idx) > 0.2*nAtoms)
%             %% too many failure likely to imply the underlying system is bad!
%             fprintf(fid,'Quit refinedistances!\n');
%             Aest = Aorg;
%             solveOK = inf;
%             break;
%         end
%     elseif (norm(Anew,'fro') > 10*norm(Astart,'fro')) | (objend > objstart)
%         Anew = Astart;
%         solveOK = 0;
%         fprintf(fid,'***** refined solution is bad, discarded\n');
%     end
%     if (solveOK==1)
%         if isempty(idxbad)
            Aest = Anew;
%         else
%             idxbad = [idxbad,idxrun(expansion_idx)];
%             Aest = zeros(size(Aorg));
%             Aest(:,setdiff(idxall,idxbad)) = Anew(:,setdiff([1:nAtoms],expansion_idx));
%             Aest(:,idxbad) = Aorg(:,idxbad);
%         end
        break;
%     else
%         beta_param = 0;
%         idxbad = [idxbad,idxrun(expansion_idx)]; %%idxbad, idxrun are for the global system
%         Aest = zeros(size(Aorg));
%         Aest(:,setdiff(idxall,idxbad)) = Anew(:,setdiff([1:nAtoms],expansion_idx));
%         Aest(:,idxbad) = Aorg(:,idxbad);
%     end
end
%%
info.objstart = objstart;
info.objend   = objend;
info.grad     = norm(Agrad,'fro')/nAtoms;
info.iter     = iter;
info.sepratio = compsepratio(Aest,D);
info.solveOK  = solveOK;
info.expansion_ratio = expansion_ratio;
fprintf(fid,'refinedistances.m: Iter = %4.0f,',iter);
fprintf(fid,' objstart = %3.2e, objend = %3.2e\n',...
    info.objstart,info.objend);
end
%%fprintf(fid,' grad = %3.2e\n',info.grad);
%%******************************************************
%% sum_{ij} (norm(xij)-dij)^2
%%******************************************************

function v = objfunc(Aest,d,eij,beta)

Aij = Aest * eij;
c   = sqrt(sum(Aij.*Aij))';
%   v   = norm((c-d))^2;
  
  weightVec=sqrt(1./d);
%   [B IX]=sort(weightVec,'descend');
%   cutoff=B(size(Aest,1));
%   weightVec(weightVec>cutoff)=cutoff;
v   = norm((c-d).*weightVec)^2;

if (abs(beta))
    v = v - beta*norm(Aest,'fro')^2;
end
end
%%******************************************************
%%******************************************************

function g = objgrad(Aest,d,eij,beta)

nEdges = length(d);
Aij    = Aest * eij;
c      = sqrt(sum(Aij.*Aij))' + eps;
%  pull   = 1 - d./c;

  weightVec=sqrt((1./d).^2);
%   [B IX]=sort(weightVec,'descend');
%   cutoff=B(size(Aest,1));
%   weightVec(weightVec>cutoff)=cutoff;
pull   = weightVec.*(1- d./c);


t = Aij * spdiags(2*pull,0,nEdges,nEdges);
g = eij*t';
g = g';
if (abs(beta))
    g = g - (2*beta)*Aest;
end
end
%%******************************************************
%%************************************************************
%% compute separation ratio
%%************************************************************
  function [sepratio,distratio] = compsepratio(Aest,D)
  
  [ii,jj,dd] = find(D); 
  nEdges = length(dd); 
  nAtoms = size(D,1); 
  eij = sparse(ii,1:nEdges,1,nAtoms,nEdges) ...
        - sparse(jj,1:nEdges,1,nAtoms,nEdges);
  Aestij = Aest*eij; 
  dest = sqrt(sum(Aestij.*Aestij)); 
  if (size(dest,1) ~= size(dd,1)); dest = dest'; end
  distratio = dest./dd; 
  sepratio = mean(distratio); 
  end
 %%************************************************************
 
