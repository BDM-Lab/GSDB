function [RHO,PVAL ] = Spearman_corr( dist1,dist2)
%% This finds the correlation between
% convert the matrix to a columns vector
% Test to observe change in perfomnace by removing  zero: knkown observed
  
    n =length(dist1);
    Dist1 = [];
    Dist2 = [];
    
    for i=1:n    
        for j=1:n    
             if(dist1(i,j)~=0)
                  Dist1 = [Dist1; dist1(i,j)]; 
                  Dist2 = [Dist2; dist2(i,j)];                 
              
             end
        end   
    end    
    
    
    % spearman correlation
 
    [RHO,PVAL] = corr(Dist1,Dist2,'Type','Spearman');
    %The spearman correlation rho
    fprintf('spearman correlation = %f\n',RHO);
  
end

