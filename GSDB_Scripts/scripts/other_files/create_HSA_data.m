%% Create a matrix data for the sparse dataset
path = '/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_normalized/GM12878/KR_1mb';

Resolution = 1000000;

 for i = 1:23
     fprintf('Processing chromosome %d............ \n', i);
     path1 = strcat(path,'/chr',int2str(i),'_matrix.txt');
     path2 = strcat(path,'/chr',int2str(i),'_HSA.txt');
     
     Data1 = dlmread(path1);
    
	 newD=[];
	 s=0; 
	for k = 1:length(Data1)	
	 e=k*Resolution;
	 D = [s e];
	 newD=[newD; D];
	 s=e;
	end         
	 Data1=[newD,Data1];   
     fprintf('file name  = %s\n', path2);	 
	 fileID = fopen(path2,'w');    
	 for j = 1:length(Data1(:,1))     
		fprintf(fileID, '%d  %d  %7.3f ',Data1(j,1:end));  
	  fprintf(fileID, '\n');   
	 end
          
     fclose(fileID);
 end
 
 disp('Completed Successfully');