function centroid = computecentroid(A)

len = size(A,2); 
centroid = (A*ones(len,1))/len; 

