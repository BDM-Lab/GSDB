%%
%% readmoive.m
%%

function [At,C,b,p,q,Gamma,Omega,samplerate,transposed] = readmovie(randstate)

if ~exist('randstate')
    randstate = 0;
end
randstate_old = rand('state');
rand('state',randstate);

M = load('u100k.data');
p=943; q=1682;
%  M = load('1Mratings.dat');
%  p = 6040; q = 3900;
ii = M(:,1);
jj = M(:,2);
rate = M(:,3);
m = length(rate);
Gamma = p*(jj-1) + ii; % here is a vector

if (p>q)
    tmp = p; p =q; q =tmp;
    ii = M(:,2);
    jj = M(:,1);
    transposed = 1;
else
    transposed = 0;
end

D = spconvert([ii,jj,rate; p,q,0]);

I = zeros(m,1);
J = zeros(m,1);
b= zeros(m,1);
count = 0;
if (transposed)
     for j=1:q
          indxrate = find(D(:,j)>0);
          indxrate = indxrate';
          for i = indxrate
              if (rand(1) < 0.5)
                  count = count+1;
                  I(count) = i;
                  J(count) = j;
                  b(count) = D(i,j);
              end
          end
     end
else
     for i=1:p
          indxrate = find(D(i,:)>0);
          % indxrate = indxrate';
          for j = indxrate
              if (rand(1) < 0.5)
                 count = count+1;
                 I(count) = i;
                 J(count) = j;
                 b(count) = D(i,j);
              end
          end
     end
end
mm = count;
I =  I(1:mm); 
J =  J(1:mm);
b = b(1:mm);


II = I + p*(J-1);   Omega = II;   % selected indices
JJ = [1:mm]';
At = spconvert([II,JJ,ones(mm,1); p*q,mm,0]);
C = zeros(p,q);

idx = setdiff(Gamma,Omega);
samplerate = D(idx);

rand('state',randstate_old);

      
