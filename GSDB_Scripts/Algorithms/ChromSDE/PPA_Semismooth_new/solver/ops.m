%%******************************************************************
%% ops:  
%%
%%   Z = ops(X,operand,Y,alpha); 
%%
%%  INPUT:        X = a matrix or a scalar
%%                    or a CELL ARRAY consisting only of matrices 
%%          operand = sym, transpose, triu, tril,
%%                    real, imag, sqrt, abs, max, min, nnz,
%%                    spdiags, ones, norm, sum, row-norm, 
%%                    rank1, rank1inv, inv
%%                    +,  -, *, .*,  ./, .^ 
%%     Y (optional) = a matrix or a scalar 
%%                    or a CELL ARRAY consisting only of matrices
%% alpha (optional) = a scalar
%%                    or the variable blk. 
%%******************************************************************

   function Z = ops(X,operand,Y,alpha); 

   spdensity = 0.5; 

   if (nargin == 2)
      if strcmp(operand,'sym'); 
         if ~iscell(X); 
            [m,n] = size(X); 
            if (m == n); 
               Z = (X+X')/2; 
            elseif (n == 1);  
               Z = X;
            else 
               error('X must be square matrix or a column vector'); 
            end 
         else
            Z = cell(size(X)); 
            for p = 1:length(X); 
                [m,n] = size(X{p}); 
                if (m == n); 
                   Z{p} = (X{p}+X{p}')/2; 
                elseif (n == 1); 
                   Z{p} = X{p};
                else; 
                   error('X{p} must be square matrix or a column vector'); 
                end 
            end
         end
      elseif strcmp(operand,'triu') | strcmp(operand,'tril'); 
         if ~iscell(X); 
            if (size(X,2) == size(X,2)); 
               eval(['Z = ',operand,'(X);']);
            elseif (size(X,2) == 1); 
               eval(['Z = X;']);
            else 
               error('X must be square matrix or a column vector'); 
            end 
         else
            Z = cell(size(X)); 
            for p = 1:length(X); 
                if (size(X{p},1) == size(X{p},2)); 
                   eval(['Z{p} = ',operand,'(X{p});']); 
                elseif (size(X{p},2) == 1); 
                   eval(['Z{p} = X{p};']);
                else; 
                   error('X{p} must be square matrix or a column vector'); 
                end 
            end
         end
      elseif strcmp(operand,'max') | strcmp(operand,'min') | ...
             strcmp(operand,'sum'); 
         if ~iscell(X); 
            eval(['Z = ',operand,'(X);']);  
         else
            Z = []; 
            for p = 1:length(X); 
                eval(['Z = [Z  ',operand,'(X{p})','  ];']);  
            end
         end 
         eval(['Z = ',operand,'(Z);']);  
      elseif strcmp(operand,'norm');
         if ~iscell(X); 
            Z = full(sqrt(sum(sum(X.*X))));
         else
            Z = 0; 
            for p = 1:length(X); Z = Z + sum(sum(X{p}.*X{p})); end;
            Z = sqrt(Z); 
         end
      elseif strcmp(operand,'trace');
         if ~iscell(X); 
            if (size(X,2)==1); Z = sum(X); else; Z = sum(diag(X)); end
         else
            Z = 0; 
            for p = 1:length(X); 
               if (size(X{p},2)==1)
                  Z = Z + sum(X{p}); 
               else
                  Z = Z + sum(diag(X{p})); 
               end
            end
         end
      elseif strcmp(operand,'inv') 
         if ~iscell(X);
            [m,n] = size(X); n2 = n*n; 
            Z = inv(X); 
            if (nnz(Z) > spdensity*n2) & issparse(Z); 
               Z = full(Z); 
            elseif (m > 1 & n == 1);
               Z = 1./X; 
               if (nnz(Z) > spdenssity*n) & issparse(Z); 
                  Z = full(Z); 
               end
            end
         else
            Z = cell(size(X)); 
            for p = 1:length(X) 
               [m,n] = size(X{p}); n2 = n*n;
               if (m == n); 
                  Z{p} = inv(X{p}); 
                  if (nnz(Z{p}) > spdensity*n2) & issparse(Z{p}); 
                     Z{p} = full(Z{p});  
                  end
               elseif (m > 1 & n == 1)
                  Z{p} = 1./X{p}; 
                  if (nnz(Z{p}) > spdensity*n) & issparse(Z{p}); 
                     Z{p} = full(Z{p}); 
                  end
               end
            end
         end
      elseif strcmp(operand,'getM')
         if ~iscell(X)
            Z = size(X,1);
         else
            for p = 1:length(X); Z(p) = size(X{p},1); end;
            Z = sum(Z); 
         end 
      elseif strcmp(operand,'nnz')
         if ~iscell(X); 
            Z = nnz(X);  
         else
            for p = 1:length(X)
               Z(p) = nnz(X{p}); 
            end;
            Z = sum(Z);  
         end 
      elseif strcmp(operand,'ones');
         if ~iscell(X)
            Z = ones(size(X));
         else
            Z = cell(size(X)); 
            for p = 1:length(X)
                Z{p} = ones(size(X{p}));
            end
         end
      elseif strcmp(operand,'zeros')
         if ~iscell(X)
            Z = sparse(size(X,1),size(X,2));
         else 
            Z = cell(size(X)); 
            for p = 1:length(X)
               Z{p} = sparse(size(X{p},1),size(X{p},2));
            end
         end
      elseif strcmp(operand,'identity');
         blk = X; 
         Z = cell(size(blk,1),1); 
         for p = 1:size(blk,1)
            pblk = blk(p,:); n = sum(pblk{2});  
            if strcmp(pblk{1},'s')
               Z{p} = speye(n,n);
            elseif strcmp(pblk{1},'q')  
               s = 1+[0, cumsum(pblk{2})];
               len = length(pblk{2});
               Z{p} = zeros(n,1);
               Z{p}(s(1:len)) = ones(len,1);
            elseif strcmp(pblk{1},'l')
               Z{p} = ones(n,1); 
            elseif strcmp(pblk{1},'u')
               Z{p} = zeros(n,1);   
            end
         end
      elseif strcmp(operand,'row-norm')
         if ~iscell(X)
            if (size(X,2) == size(X,1)); 
               Z = sqrt(sum((X.*conj(X))'))';
            elseif (size(X,2) == 1);  
               Z = abs(X); 
            end;
         else
            Z = cell(size(X)); 
            for p = 1:length(X)
               if (size(X{p},2) == size(X{p},1))
                  Z{p} = sqrt(sum((X{p}.*conj(X{p}))'))'; 
               elseif (size(X{p},2) == 1)  
                  Z{p} = abs(X{p}); 
               end
            end
         end       
      elseif strcmp(operand,'sqrt') | strcmp(operand,'abs') | ...
         strcmp(operand,'real') | strcmp(operand,'imag') | ...
         strcmp(operand,'transpose') | strcmp(operand,'log');
         if ~iscell(X); 
            eval(['Z = ',operand,'(X);']);  
         else
            Z = cell(size(X)); 
            for p = 1:length(X); 
               eval(['Z{p} = ',operand,'(X{p});']);  
            end
         end
      end
   end
%%
   if (nargin == 3) 
      if strcmp(operand,'spdiags');
         if ~iscell(Y)
            [m,n] = size(Y); 
            if (m == n); 
               Z = spdiags(X,0,m,n);
            else
               Z = X;
            end;
         else 
            Z = cell(size(Y)); 
            for p = 1:length(Y);
               [m,n] = size(Y{p}); 
               if (m == n); 
                  Z{p} = spdiags(X{p},0,m,n);
               else;
                  Z{p} = X{p};
               end
            end
         end
      elseif strcmp(operand,'inprod')
         if ~iscell(X) & ~iscell(Y)
    	    Z = (Y'*X)';
	 elseif iscell(X) & iscell(Y)
    	    Z = zeros(size(X{1},2),1); 
   	    for p=1:length(X)
	        Z = Z + (Y{p}'*X{p})'; 
            end
         end
      elseif strcmp(operand,'+') | strcmp(operand,'-') | ...
             strcmp(operand,'/') | strcmp(operand,'./') | ...
             strcmp(operand,'*') | strcmp(operand,'.*') | ...
             strcmp(operand,'.^')
         if (~iscell(X) & ~iscell(Y)) 
            eval(['Z = X',operand,'Y;']); 
         elseif (iscell(X) & iscell(Y))
            Z = cell(size(X)); 
            for p = 1:length(X) 
 	           if (size(X{p},2) == 1) & (size(Y{p},2) == 1) & ... 
                  (strcmp(operand,'*') | strcmp(operand,'/')) 
                  eval(['Z{p} = X{p}.',operand,'Y{p};']); 
               else
                  eval(['Z{p} = X{p} ',operand,'Y{p};']);                   
               end;
            end; 
         elseif (iscell(X) & ~iscell(Y))
            Z = cell(size(X)); 
            for p = 1:length(X); 
                eval(['Z{p} = X{p}',operand,'Y;']); 
            end
         elseif (~iscell(X) & iscell(Y))
            Z = cell(size(Y)); 
            for p = 1:length(Y); 
                eval(['Z{p} = X',operand,'Y{p};']); 
            end
         end
      else
         error([operand,' is not available, check input arguments']); 
      end
   end
%%
   if (nargin == 4) 
      if strcmp(operand,'+') | strcmp(operand,'-')
         if ~iscell(X) & ~iscell(Y)
            eval(['Z = X',operand,'alpha*Y;']); 
         elseif (iscell(X) & iscell(Y))
            Z = cell(size(X)); 
            for p = 1:length(X) 
               eval(['Z{p} = X{p}',operand,'alpha*Y{p};']); 
            end
         else
            error('X, Y are different objects'); 
         end
      else
         error([operand,' is not available']); 
      end
   end
%%******************************************************************
