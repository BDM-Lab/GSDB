%%
%% assume p<= q
%%

function [At,C,b,M,Omega] = randMC(p,q,r,ratio,noiseLevel,randstate)

if ~exist('noiseLevel'); noiseLevel = 0; end 
if ~exist('randstate'); randstate = 0; end 

  randnstate_old = randn('state');
  randstate_old = rand('twister');
 
  randn('state',randstate);
  rand('twister',randstate);

  ML = randn(p,r); 
  MR = randn(q,r); 
  M = ML*MR'; 
%%
  dr = r*(p+q-r);
  ntotal = p*q; 
  m  = ratio*dr;
  Omega = zeros(m,1);  count = 0;
  prob = min(1,m/ntotal);
  for j=1:q
      tmp = rand(p,1);
      idx = find(tmp<prob);
      len = length(idx);
      Omega(count+[1:len]) = idx + p*(j-1);
      count = count + len;
  end
  mm = count;
  Omega = Omega(1:mm);
  vv = M(Omega);
  JJ = [1:mm]';
  %%
  At = spconvert([Omega,JJ,ones(mm,1); ntotal,mm,0]);
  if (noiseLevel > 0)
     randvec = randn(mm,1); 
     b = vv + noiseLevel*norm(vv)*randvec/norm(randvec);  
  else
     b = vv; 
  end
  C = zeros(p,q);
  
  randn('state',randnstate_old);
  rand('twister',randstate_old);
  
  
