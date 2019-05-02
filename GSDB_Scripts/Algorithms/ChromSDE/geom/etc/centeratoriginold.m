function A = centeratorigin(A)

npts = size(A,2); 
center = A*ones(npts,1)/npts; 
A = A - center*ones(1,npts); 

