function [binAnno,FreqMat,normFreqMat]=readpipeline_output(prefix)
spMat= importdata([prefix '.n_contact'], '\t', 1);
sel=spMat.data(:,4)>0;%need more than 2 pet count support
FreqMat=spconvert(spMat.data(sel,[1 2 4]));
FreqMat=max(FreqMat,FreqMat'); %
normF=spMat.data(sel,4)./spMat.data(sel,3); 
normFreqMat=spconvert([spMat.data(sel,[1 2]),normF]);
normFreqMat=max(normFreqMat,normFreqMat');
binData=importdata([prefix '.cbins'], '\t', 1);
binAnno=binData.data(:,[2 3 4 1]);
binAnno=checkMissingBin(binAnno);
end


function binAnno=checkMissingBin(binAnno)

for m=1:size(binAnno,1)
    if(m~=binAnno(m,4))
        binAnno=[binAnno(1:m-1,:);[0 0 0 m];binAnno(m:end,:)];
    end
    
end


end