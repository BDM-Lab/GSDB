function A = translate(A,v)

[nDimensions,nAtoms] = size(A);
A = A + v * ones(1,nAtoms);

