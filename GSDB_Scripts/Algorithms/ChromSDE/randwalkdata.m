
function [binAnno FreqMat XX]=randwalkdata(n,noiserate)
rand('twister',1234);
addpath('./helperfunctions');
XX=randwalk(n);
plot3(XX(1,:),XX(2,:),XX(3,:));
hold on
plot3(XX(1,:),XX(2,:),XX(3,:),'r.');
hold off
drawnow

FreqMat=points2FreqMat(XX,noiserate);
binAnno=[ones(size(FreqMat,1),1),(1:size(FreqMat,1))'];
end


function FreqMat=points2FreqMat(X,noiselevel)

DD=sqrt(dissimilar(X,inf,50));
[a b c]=find(triu(DD));
FreqMat=zeros(size(DD));
for i=1:length(a)
   lambda=1/c(i); 
   lambda=lambda*(1+(rand-0.5)*noiselevel*2);
   F=lambda;
   FreqMat(a(i),b(i))=F;
end
sumcount=sum(FreqMat(:));
FreqMat=FreqMat/sumcount*1000000;
FreqMat=max(FreqMat,FreqMat');
end

function Y=randwalk(n)
npoints =n;
dt = 1;
bm = cumsum([zeros(1, 3); dt^0.5*randn(npoints-1, 3)]);

Y=bm';

end