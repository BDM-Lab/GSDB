function [binAnno FreqMat XX]=randdata(n,noiserate)
rand('twister',1234);
XX=rand(3,n);
addpath('./helperfunctions');
plot3(XX(1,:),XX(2,:),XX(3,:));
hold on
plot3(XX(1,:),XX(2,:),XX(3,:),'r.');
hold off
drawnow

FreqMat=points2FreqMat(XX,noiserate);
binAnno=[ones(size(FreqMat,1),1),(1:size(FreqMat,1))'];
end

function FreqMat=points2FreqMat(X,noiselevel)
nn=50;

DD=sqrt(dissimilar(X,inf,round(nn)));
[a b c]=find(triu(DD));
FreqMat=zeros(size(DD));
for i=1:length(a)
   lambda=(1/c(i)); 
   lambda=lambda*(1+(rand-0.5)*noiselevel*2);
   F=lambda;
   FreqMat(a(i),b(i))=F;
end
sumcount=sum(FreqMat(:));
FreqMat=FreqMat/sumcount*1000000;
FreqMat=max(FreqMat,FreqMat');
end

