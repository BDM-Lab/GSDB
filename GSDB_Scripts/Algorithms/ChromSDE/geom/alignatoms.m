%%**********************************************************************
%% alignatoms Aligns the estimated atom positions to the true positions.
%%
%% Info's fields
%% -------------
%% rmsd
%% Transformation : transform(Aest,Transformation) aligns Aest to A.
%% Aest
%%
%% Info = alignatoms(A,Aest)
%%
%% Created       :  6 Jul 07
%% Last modified :  6 Jun 08
%%**********************************************************************

  function Info = alignatoms(A,Aest)
   npts = size(Aest,2); 
    centroidB     = computecentroid(Aest);
     centroidA      = computecentroid(A);
     Aest              = Aest - centroidB*ones(1,npts); 
     
     A              = A - centroidA*ones(1,npts); 
  %scale first
  scale=norm(A)/norm(Aest);
  Aest=scale*Aest;
  options = 2; %% better 
  if (options==1)
     [dummy,B,Transformation] = procrustes(A',Aest'); 
     B = B'; 
     Info.Aest = B;
     Info.rmsd = computermsd(A,B,false);
     Transformation.procrustes = 1; 
     Transformation.T    = Transformation.T'; 
     Transformation.c    = Transformation.c(1,:)'; 
     Info.Transformation = Transformation; 
  else
      
     B = Aest; 
     
  

     rotationMatrix = orthogonalprocrustes(A,B);
     B              = rotationMatrix *B; 
     
     Info.rmsd                      = computermsd(A,B,false);
     Info.Transformation.procrustes = 0; 
     Info.Transformation.aTranslate = -centroidB;
     Info.Transformation.bRotate    = rotationMatrix*scale;
     Info.Transformation.cTranslate = centroidA;
     Info.Aest                      = B + centroidA*ones(1,npts); 
     
  end
%%**********************************************************************
%%**********************************************************************
%% orthogonalprocurstes Solve the orthogonal Procrustes problem.
%%
%% Solves the orthogonal Procrustes problem:
%%            Find U such that U^TU = I, and U minimizes ||A-UB||_F.
%%
%% U = orthogonalprocrustes(A,B)
%%**********************************************************************

  function U = orthogonalprocrustes(A,B)

  [W,S,Z] = svd(A*B');
  U = W * Z';
%%**********************************************************************
%%**********************************************************************
  function centroid = computecentroid(A);

  len = size(A,2);
  centroid = (A*ones(len,1))/len;
%%**********************************************************************
