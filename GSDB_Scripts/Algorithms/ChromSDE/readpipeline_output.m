function [binAnno,normFreqMat]=readpipeline_output(prefix)
spMat= importdata([prefix '.n_contact'], '\t');
normFreqMat=spconvert(spMat);
normFreqMat=max(normFreqMat,normFreqMat');
% normF=spMat.data(sel,4)./spMat.data(sel,3); %//replace with normalized frequency
% normFreqMat=spconvert([spMat.data(sel,[1 2]),normF]);
% normFreqMat=max(normFreqMat,normFreqMat');
%binData=importdata([prefix '.cbins'], '\t');
%binAnno=binData.data(:,[2 3 4 1]);
binAnno=importdata([prefix '.cbins'], '\t'); %Already in the format [2 3 4 1] from file
binAnno=checkMissingBin(binAnno);
end


function binAnno=checkMissingBin(binAnno)

for m=1:size(binAnno,1)
    if(m~=binAnno(m,4))
        binAnno=[binAnno(1:m-1,:);[0 0 0 m];binAnno(m:end,:)];
    end
    
end


end