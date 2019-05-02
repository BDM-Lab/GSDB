%%******************************************************
%% transform Apply a transformation.
%%
%% points = transform(points,Transformation)
%%******************************************************

  function A = transform(A,Transformation)

  npoints = size(A,2);
  A = A +  Transformation.aTranslate * ones(1,npoints);
  A = Transformation.bRotate * A; 
  A = A + Transformation.cTranslate * ones(1,npoints);
%%******************************************************
