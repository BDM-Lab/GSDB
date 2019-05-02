function [error_term Pest2 ConsensusIndex]=ChromSDE_knownAlpha(heatmapMatrix,alpha,method_type)
         D_raw=zeros(size(heatmapMatrix));
         D_raw(heatmapMatrix>0)=(heatmapMatrix(heatmapMatrix>0).^(-alpha));


        scale_d=max(max(D_raw(D_raw>0)));
        
         D_raw=D_raw/scale_d;
        D=D_raw;

%%start config optimization
 rhofactor = 0.01; 

if method_type
        addpath('./PPA_Semismooth_new/solver');
        addpath('./PPA_Semismooth_new/util');
         addpath('./PPA_Semismooth_new/runexpt');
         addpath('./PPA_Semismooth_new/runexpt/PDB');
         addpath('./PPA_Semismooth_new/solver/mexfun');
         DD = D.*D; n = length(D);
         
         Rank = 3;
         DD = triu(DD,1)+triu(DD,1)';
         H = spones(DD); % equal weight    
         [ai bi ci]=find(triu(DD));
         
 H_options = 1; %% 1 is better
         if (H_options==1) 
             weightVec=sqrt(1./ci);
              [B IX]=sort(weightVec,'descend');
              cutoff=B(size(DD,1));
              weightVec(weightVec>cutoff)=cutoff;
%               weightVec=ones(size(weightVec));
            H = spconvert([ai,bi,weightVec;n,n,0]);
            H = H+H';  
         end
         [blk,At,b,numeq,PrecondData] = RKESDPdata(DD,H);
         eigsopts.disp  = 0;
         eigsopts.tol   = 1e-6;
         eigsopts.maxit = 200;
         Atb = Atyfunsym(blk,At,b,1);
         [dummy,Lam,flag] = eigs(Atb{1},3,'LM',eigsopts);
         rhomax = max(abs(diag(Lam)));
         rho = -rhofactor*rhomax*nnz(DD)/(size(DD,1)^2);
         OPTIONS.tol  = 1e-6; 
         OPTIONS.sig0 = 100; 
         OPTIONS.continuation  = 1; 
         OPTIONS.precond       = 3; 
         OPTIONS.Newton_alteroption = 0;
         OPTIONS.maxitpsqmr    = 600;  % PSQMR is very expensive
         OPTIONS.plotyes = 0;
         OPTIONS.printlevel=0;
         if (OPTIONS.precond == 3) 
            OPTIONS.PrecondData = PrecondData;
         end
        
 C_options=1; %% 1 is better
         if (C_options==0)
            C{1} = sparse(n,n); 
         elseif (C_options==1)
            cnt = 0; 
            row=zeros(length(ai)*4,1);
            col=zeros(length(ai)*4,1);
            vv=zeros(length(ai)*4,1);
           normci=ci/mean(ci);
            for kk = 1:length(ai)          
                  row(cnt+[1:4],1) = [ai(kk),ai(kk),bi(kk),bi(kk)];
                  col(cnt+[1:4],1) = [ai(kk),bi(kk),ai(kk),bi(kk)];
                   weight= abs(rho)/size(DD,1);
                  vv(cnt+[1:4],1) = weight*[1,-1,-1,1];
                  cnt = cnt+4;           
            end   
            C{1} = spconvert([row,col,vv;n,n,0]);
         end
         [obj,XX,yy,ZZ,PPAInfo,runhist] = PPASolversym(blk,At,b,numeq,C,rho,OPTIONS);
     
         Record.rhofactor = rhofactor;
         Record.rho = rho;
         Record.constraints = length(bi);
         Record.cputime = runhist.cputime(end);
         X0 = XX{1}; 
         e = ones(n,1);

         [V, eigX] = mexeig(full(X0));
         if (size(eigX,1)==size(eigX,2)); eigX = diag(eigX); end
                  %%%%%%%satisfyRatio%%%%%%%%%%
         K=full(X0);
          D = triu(D,1)+triu(D,1)';
         [a,b,c] = find(D);
         sizeD=size(D);
           Iaa=sub2ind(sizeD,a,a);
         Iba=sub2ind(sizeD,b,a);
        Iab=sub2ind(sizeD,a,b);
          Ibb=sub2ind(sizeD,b,b);
       % satisfiedRatio=sum(abs(D(Iab)-sqrt((K(Iaa)-K(Iab)-K(Iba)+K(Ibb))))<0.2*D(Iab))/length(Iba);
       satisfiedRatio=sum(min(D(Iab)./sqrt((K(Iaa)-K(Iab)-K(Iba)+K(Ibb))),sqrt((K(Iaa)-K(Iab)-K(Iba)+K(Ibb)))./D(Iab)))/length(Iba);
          %%%%%%%satisfyRatio%%%%%%%%%%
         [dummy,idx] = sort(eigX,'descend');
         eigX  = full(eigX(idx)); V = V(:,idx);
          ConsensusIndex=sum(((eigX(1:Rank))))/sum(((eigX(eigX>0))))*satisfiedRatio;
         sprintf('Consensus Ratio: %f',ConsensusIndex)
         Pest0 = diag(sqrt(eigX(1:Rank)))*V(:,1:Rank)'; 
         posMatrix=[];
         

else
        rmpath('./PPA_Semismooth_new/solver');
        
   addpath(genpath('./SDPT3-4.0/')); 
 addpath(genpath('./yalmip'));


 [Pest0,ConsensusIndex]=linear_SDP(D,rhofactor);
  Pest0=Pest0';
       

end

        [Pest2,Info] = refinepositions(Pest0,D,0.0001);
          
        
error_term=FitFreq(squareform(pdist(Pest2')),heatmapMatrix,-alpha,0.01);
sprintf('alpha:%f  error:%f',alpha,error_term)

end