
%%*******************************************************************
%% [blk,At,b,M,Omega] =  randMCsym(p,r,ratio,noiseLevel,randstate);
%% ratio = m/dr;
%%*******************************************************************

function [blk,At,b,M,Omega] = randMCsym(p,r,ratio,noiseLevel,randstate)

if ~exist('noiseLevel'); noiseLevel = 0; end
if ~exist('randstate'); randstate = 0; end

randnstate_old = randn('state');
randstate_old = rand('state');

randn('state',randstate);
rand('state',randstate);


ML = randn(p,r);
%% normML = sqrt(sum(ML.*ML))';
%% put the points on the unit sphere.
%% ML = ML*spdiags(1./normML,0,p,p);

M = ML*ML';
M = 0.5*(M+M');
%%
dr = r*p - (r-1)*r/2;
ntotal = p*(p+1)/2;
m  = ratio*dr;
ii = zeros(m,1);
jj = zeros(m,1);
% vv = zeros(m,1);
prob = min(1,m/ntotal);
%%
count = 0;
for j = 1:p
    tmp = rand(j,1);
    indx = find(tmp<prob);
    len = length(indx);
    ii(count+[1:len]) = indx;
    jj(count+[1:len]) = j*ones(len,1);
    count = count + len;
    %     for i = 1:j
    %         if (rand(1) < prob)
    %             count = count+1;
    %             ii(count) = i;
    %             jj(count) = j;
    %             vv(count) = M(i,j);
    %         end
    %     end
end
mm = count;
ii = ii(1:mm);
jj = jj(1:mm);
indx2 = p*(jj-1) + ii;
% vv = vv(1:mm);
vv = M(indx2);
Omega = spconvert([ii,jj,ones(mm,1); p,p,0]);
Omega = spones(Omega+Omega');
%%
blk{1,1} = 's'; blk{1,2} = p;
options = 1;
if (options==1)
    II = ii + (jj-1).*jj/2;
    JJ = [1:mm]';
    VV = (1/sqrt(2))*ones(mm,1);
    idx = find(ii==jj);
    VV(idx) = ones(length(idx),1);
    At{1} = spconvert([II,JJ,VV; p*(p+1)/2,mm,0]);
elseif (options==2)
    Acell = cell(1,mm);
    for k = 1:mm
        ik = ii(k); jk = jj(k);
        if (ik < jk);
            Acell{k} = spconvert([ik,jk,0.5; jk,ik,0.5; p,p,0]);
        elseif (ik==jk)
            Acell{k} = spconvert([ik,jk,1; p,p,0]);
        end
    end
    At = svec(blk,Acell,1);
end
if (noiseLevel > 0)
    randvec = randn(mm,1);
    b = vv + noiseLevel*norm(vv)*randvec/norm(randvec);
else
    b = vv;
end
%%
randn('state',randnstate_old);
rand('state',randstate_old);
%%*******************************************************************
