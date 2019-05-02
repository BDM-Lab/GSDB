function [posMatrix,ConsensusIndex,alpha,runningTime]=ChromSDE_single(binAnno, heatmapMatrix, method_type)

tic;
mappableRegion=sum(heatmapMatrix>0)>2; %at least connect to two points
selchrom=(binAnno(:,1)<100)&mappableRegion';

binAnno=binAnno(selchrom,:);
heatmapMatrix=heatmapMatrix(selchrom,:);
heatmapMatrix=heatmapMatrix(:,selchrom);

  
FreqMat=heatmapMatrix;
binsize=abs(binAnno(1,2)-binAnno(2,2))

%[num_parts,comps]=graphconncomp(sparse(FreqMat),'Weak', 1);
[num_parts,comps] = conncomp(FreqMat);
%%make sure connected   
if num_parts>1
    minbinFreq=max(max(FreqMat(FreqMat>0)));
       for i=2:size(binAnno,1)
            if binAnno(i,2)-binAnno(i-1,2)==binsize
                if FreqMat(i-1,i)< minbinFreq
                    minbinFreq=FreqMat(i-1,i);
                end
            end
        end
        
        for i=2:size(binAnno,1)
            if  FreqMat(i-1,i)==0 &&  binAnno(i,1)==binAnno(i-1,1)&& comps(i-1)~=comps(i)
                   FreqMat(i-1,i)=minbinFreq/(abs(binAnno(i,2)-binAnno(i-1,2))/binsize);
                   FreqMat(i,i-1)=FreqMat(i-1,i);
            end
        end
end


logalpha= fminbnd(@(x) ChromSDE_knownAlpha(FreqMat,exp(x),method_type),log(0.1),log(3),optimset('TolX',0.001,'Display','off','MaxIter',10));
alpha=exp(logalpha);
[error_term Pest2 ConsensusIndex]=ChromSDE_knownAlpha(FreqMat,alpha,method_type);
 
    runningTime = toc;
     
         for pass = [1:2]
            if (pass==1); continue;  end
            genenumber = binAnno(:,1);
            idxstart = 1+[0; find(diff(genenumber)); length(genenumber)];
            for kk = 1:length(idxstart)-1
               randcolor = rand(1,3); 
               for k = idxstart(kk):idxstart(kk+1)-2
                  if (pass==1)
                     x = Pest0(:,k); y = Pest0(:,k+1);                       
                  else
                     x = Pest2(:,k); y = Pest2(:,k+1);                                                
                  end
                  % h = plot3([x(1); y(1)],[x(2); y(2)],...
                  %  [x(3); y(3)],'linewidth',3,'color',randcolor);
                  % hold on;
                  % plot3([x(1); y(1)],[x(2); y(2)],[x(3); y(3)],'.k','markersize',2);
               end    
            end 
            % grid; axis('square'); 
            % hold off;
         end
         

posMatrix=[-1 binAnno(1,1) 1 alpha];     
bigend=10000000000;
for i=1:size(Pest2,2)
    posMatrix= [posMatrix;binAnno(i,1)*bigend+binAnno(i,2),Pest2(1,i),Pest2(2,i),Pest2(3,i)];    
end

%%*******************************************************************
end
