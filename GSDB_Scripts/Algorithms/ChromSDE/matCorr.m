function [spear,pear,rmse]=matCorr(M1,M2)


if size(M1,1)~=size(M1,2)
     if size(M1,1)>size(M1,2)
         M1=squareform(pdist(M1));
         
     else
         M1=squareform(pdist(M1'));
     end
end



if size(M2,1)~=size(M2,2)
     if size(M2,1)>size(M2,2)
         M2=squareform(pdist(M2));
         
     else
         M2=squareform(pdist(M2'));
     end
end

%ignore the adjection bins and diagnal bins, second adjection bins
ii=sub2ind(size(M1),1:size(M1,1),1:size(M1,1));
ij=sub2ind(size(M1),1:size(M1,1)-1,2:size(M1,1));
ji=sub2ind(size(M1),2:size(M1,1),1:size(M1,1)-1);
ij2=sub2ind(size(M1),1:size(M1,1)-2,3:size(M1,1));
ji2=sub2ind(size(M1),3:size(M1,1),1:size(M1,1)-2);
M1(ii)=0;
M1(ij)=0;
M1(ji)=0;
M1(ij2)=0;
M1(ji2)=0;
X=M1(M1~=0&M2~=0);
Y=M2(M1~=0&M2~=0);

if isempty(X) || isempty(Y)
    spear=0;
    pear=0;
    rmse=0;
else
       
spear=corr(X(:), Y(:), 'type', 'Spearman');

pear=corr(X(:), Y(:));

% First calculate the "error".
 err = (X(:)- Y(:));
 squareError = err(:).^2;
 meanSquareError = mean(squareError);
rmse = sqrt(meanSquareError);
 


end
end