
% ShRec3D script
disp (['Running Job for filename =', filename]);
Data2  = dlmread(filename);
Data=1./(Data2);
n = size(Data, 1);	
for k=1:n
	i2k = repmat(Data(:,k), 1, n);
	k2j = repmat(Data(k,:), n, 1);
	Data = min(Data, i2k+k2j);
end   
n=length(Data);
center=zeros(n,1);
for i=1:n
	for j=1:n
	center(i)=center(i)+Data(i,j)^2-1/n*(Data(j,j:n)*Data(j:n,j));
	end
end
center(center<0)=0;
center=sqrt(center/n);
distmat=1/2*(repmat(center.*center,1,n)+repmat((center.*center)',n,1)-Data.*Data);
[V,D]=eigs(distmat,4);
lambda=diag(D);
XYZ=[V(:,1)*sqrt(lambda(1)) V(:,2)*sqrt(lambda(2)) V(:,3)*sqrt(lambda(3))];
scale=100/max(max(XYZ));
XYZ=XYZ*scale;

%output pdb and image  
output_structure;
% correlate structure
[dist] = StructureDistance(XYZ,n); % distance from structure
%convert IF to distance
for j = 1:n
	for k = 1:n
		if(Data2(j,k)~=0)
			Data(j,k)=1./(Data2(j,k));
		 end
	end
end

[RHO,RHO2] = Spearman_corr(Data,dist,Data2,n); 
fileID=fopen(logout,'w');
fprintf(fileID, 'Spearman correlation Expected Dist vs. Reconstructed Dist: %f\n',RHO);
%fprintf(fileID, 'Pearson correlation Expected Dist vs. Reconstructed Dist: %f\n',RHO2);
fclose(fileID);

