function [binAnno FreqMat X]=helixdata(nn,noiselevel)
addpath('./helperfunctions');
%%helix data
step=pi/nn*10;
t = step:step:10*pi;
X=[sin(t);cos(t);t/10];
plot3(X(1,:),X(2,:),X(3,:));
hold on
plot3(X(1,:),X(2,:),X(3,:),'r.');
hold off
drawnow
DD=sqrt(dissimilar(X,inf,50)); %fix to 50 nn
[a b c]=find(triu(DD));
FreqMat=zeros(size(DD));
binAnno=[ones(size(DD,1),1),(1:size(DD,1))'];
scale=2*max(c);
for i=1:length(a)
   lambda=(scale/c(i)); 
   lambda=lambda*(1+(rand-0.5)*noiselevel*2);
   %F= poissrnd(lambda);
   F=lambda;
   FreqMat(a(i),b(i))=F;
end
sumcount=sum(FreqMat(:));
FreqMat=FreqMat/sumcount*1000000;
FreqMat=max(FreqMat,FreqMat');
end