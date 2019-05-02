%% Create a matrix data for the sparse dataset
path = '/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_normalized/GM12878/KR_250kb';

Resolution = 250000;


 for i =1:23
     fprintf('Processing chromosome %d............ \n', i);
     path1 = strcat(path,'/chr',int2str(i),'_matrix.txt');

     path3 = strcat(path,'/chr',int2str(i),'_miniMDS.bed');
     Data = dlmread(path1);	 
	
	 newD=[];
	 for j=1:length(Data)
	       if (sum(Data(j,:))==0)
		  %  fprintf('Skipped index = %d \n', j);
			continue;
		   end
		 s=(j-1)*Resolution;
		 for k=j:length(Data)
		    if (Data(j,k)==0)
			 continue;
			end
		     e=(k-1)*Resolution;
			D = [s e Data(j,k)];
			newD=[newD; D];
		end
	end
	 
	 
	 

     fprintf('file name  = %s\n', path3);		 

	 fileID1 = fopen(path3,'w');  
     for k=1: length(newD) 
          chr=['chr',int2str(i)];
          fprintf(fileID1,'%6s %12d %12d %6s %12d %12d %9.3f \n',chr,newD(k,1), newD(k,1)+Resolution, chr, newD(k,2), newD(k,2)+Resolution, newD(k,3));  		  
     end 
	 
	 fclose(fileID1);
 end
   
    
	   
 disp('Completed Successfully');