%%*************************************************************************
%% generate SDP data corresponding to
%%
%%   min_{X psd} sum_{ij} w(i,j)(<Aij,X> - dij)^2 + lam*Tr(X)
%%
%% input: DD = (npts)x(npts) dis-similar matrix
%%*************************************************************************

function  [blk,At,bb,numeq,PrecondData] = RKESDPdata(DD,WW)

  if (nargin < 2); WW = spones(DD); end
  [ii,jj,vv] = find(WW);
  if (min(vv) < 0); error('nonzero element of WW must be positive'); end
  npts = size(DD,2);
  if (size(DD,1) ~= size(DD,2)); error(' DD must be square'); end
  DD = triu(DD,1);
  DD = DD.*spones(WW);
  NZ = nnz(DD);
%%
  sqrt2 = sqrt(2);
  bb = zeros(NZ+1,1);
  ww = zeros(NZ,1);
  indx_hbar = zeros(NZ,1);
%%
  row_indx = zeros(2*NZ,1);
  col_indx = zeros(2*NZ,1);
  sqrtw = zeros(2*NZ,1);
  indx_h = zeros(2*NZ,1);
%%
  II = zeros(3*NZ,1);
  JJ = zeros(3*NZ,1);
  VV = zeros(3*NZ,1);
  count = 0; count2 = 0;
  mm = 0;
for k = 1:npts
    for j = 1:k-1
        rr = DD(j,k);
        if (rr) 
            mm = mm+1;
            sqrtwjk = sqrt(WW(j,k));
            bb(mm) = sqrtwjk*rr;
            ww(mm) = WW(j,k);
            indx_hbar(mm) = k*(k-1)/2+j; % jk
            row_indx(count2 + [1:2]) = [j; k];
            col_indx(count2 + [1:2]) = mm*ones(2,1); % for Jtranspose
            sqrtw(count2 + [1:2]) = [sqrtwjk; sqrtwjk];
            indx_h(count2 + [1:2]) = [j*(j+1)/2; k*(k+1)/2]; % [jbar; kbar];
            %%
            %% hh = h(indx_h); 
            %% J = sparse(NZ,npts); J(indx_J) = sqrtw.*sqrt(hh);
            %%
            II(count+[1:3]) = [j*(j+1)/2;(k-1)*k/2+j; k*(k+1)/2]; % [jbar; jk; kbar];
            JJ(count+[1:3]) = mm*ones(3,1);
            VV(count+[1:3]) = sqrtwjk*[1;-sqrt2;1];
            count = count + 3;
            count2 = count2 + 2;
        end
    end
end
%%
%% set e'*Y*e = 0 to center points around origin
%%
  blk{1,1} = 's'; blk{1,2} = npts;
  e = sqrt(2/npts)*ones(npts,1);
  Alast = svec(blk,e*e',1);
  At{1} = [spconvert([II,JJ,VV;npts*(npts+1)/2,mm,0]), Alast];
  bb(mm+1) = 0;
  numeq = 1;
%%
  PrecondData.ww = ww;
  PrecondData.indx_hbar = indx_hbar;
  PrecondData.row_indx = row_indx;
  PrecondData.col_indx = col_indx;
  PrecondData.sqrtw = sqrtw;
  PrecondData.indx_h = indx_h;
%************************************************************************
