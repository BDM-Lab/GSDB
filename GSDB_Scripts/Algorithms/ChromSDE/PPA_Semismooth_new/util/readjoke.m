%%
%% readjokens.m
%% read Jester Joke data set for nuclear norm minimization
%%

function [At,C,b,p,q,Gamma,Omega,samplerate,transposed] = readjoke(Data, randstate)

if ~exist('randstate')
    randstate = 0;
end
randstate_old = rand('twister');
rand('twister',randstate);

M = Data(:,[2:101]);
[p,q]= size(M);
p0 = p; M0 = M;

[ii,jj, rate]=find(M~=99);
m = length(rate);
Gamma = p*(jj-1) + ii; % here Omega is a vector

%%
if (p>q)
    tmp = p; p =q; q =tmp;
    M = M';
    transposed = 1;
else
    transposed = 0;
end

%% randomly select 10 ratings per user
I = zeros(m,1);
J = zeros(m,1);
b= zeros(m,1);
count = 0;

for j=1:q  % Since in the datset p>q, transposed = 1;
    jokenum = Data(j,1);
    prob = 10/jokenum;
    indxrate = find(M(:,j)~=99); % here indxrate is column vector
    indxrate = indxrate';
    for i = indxrate
        if (rand(1) < prob)
            count = count+1;
            I(count) = i;
            J(count) = j;
            b(count) = M(i,j);
        end
    end
end

mm = count;
I =  I(1:mm);
J =  J(1:mm);
b = b(1:mm);

II = I + p*(J-1);
JJ = [1:mm]';
At = spconvert([II,JJ,ones(mm,1); p*q,mm,0]);
C = zeros(p,q);


III = J;
JJJ = I;
Omega = III + p0*(JJJ-1);
idx = setdiff(Gamma,Omega);
samplerate = M0(idx);

rand('twister',randstate_old);








