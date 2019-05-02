%method_type=1 Quadratic SDP
%method_type=0 Linear SDP
function ChromSDE(trainBin,trainFreq,method_type)
bigend=10000000000;
figure
record=[];
posMatrix_combined=[];
chromNum=min(24,max(trainBin(:,1)));
rowNum=ceil(chromNum/5);
colNum=min(5,chromNum);
for chrom=1:chromNum
    selchrom=(trainBin(:,1)==chrom);
    if sum(selchrom)==0
        continue
    end
    binAnno=trainBin(selchrom,:);
    normFreqMat=trainFreq(:,selchrom);
    normFreqMat=normFreqMat(selchrom,:);
    subplot(rowNum,colNum,chrom);
    sprintf('chrom %f  npt:%f',chrom,size(binAnno,1))
  [posMatrix,consensusIndex,alpha]=ChromSDE_single(binAnno, normFreqMat,method_type);

    YY=posMatrix(posMatrix(:,1)>0,2:4);
    PosList=posMatrix(posMatrix(:,1)>0,1);
    D=squareform(pdist(YY));
    D_Freq=zeros(size(D));

    D_Freq(D>0)=(D(D>0).^(-1/alpha));
    [~,ia,ib] = intersect(PosList-chrom*bigend,binAnno(:,2));
    [S1,P1]=matCorr(D_Freq(ia,ia),normFreqMat(ib,ib));
     posMatrix_combined=[posMatrix_combined;posMatrix];
         title([chrom S1]);
    end
    record=[record;chrom P1 S1 consensusIndex alpha];
    

name=inputname(2);
dlmwrite([name  '.stat.txt'],sprintf('chrom\tPearson\tSpearman\tConsensusIndex\tAlpha'));
 dlmwrite([name  '.stat.txt'], record,'-append', 'precision', '%5.5g', 'delimiter', '\t');
 dlmwrite([name  '.pos'], full(posMatrix_combined), 'precision', '%10.10g');
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
 print(gcf, '-dpng', [name '.png']);
saveas(gcf,name,'fig')
pos2pdb([name  '.pos']);
end