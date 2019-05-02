function [RHO,RHO2 ] = Spearman_corr( dist1,dist2,contact,n)
%% This finds the correlation between
% convert the matrix to a columns vector
% Test to observe change in perfomnace by removing  zero: knkown observed
  
    Dist1 = [];
    Dist2 = [];
    
    for i=1:n    
        for j=1:n    
             if(contact(i,j)~=0)
                  Dist1 = [Dist1; dist1(i,j)]; 
                  Dist2 = [Dist2; dist2(i,j)];                 
              
             end
        end   
    end    
    
    
    % spearman correlation 
    [RHO,PVAL1] = corr(Dist1,Dist2,'Type','Spearman');
    %The spearman correlation rho
	fprintf('spearman correlation = %f\n',RHO);
   % Pearson correlation 
    [RHO2,PVAL2] = corr(Dist1,Dist2);
    %The Pearson correlation rho2
	%fprintf('pearson correlation = %f\n',RHO2);
end

