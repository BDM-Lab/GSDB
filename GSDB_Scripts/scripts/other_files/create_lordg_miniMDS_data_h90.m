%% Create a matrix data for the sparse dataset
path = '/storage/htc/bdm/tosin/';
sqpath= [path,'GSDB/Data/BB8015WF/hIMR90'];

fullPath = fullfile(sqpath, 'List_data');
if ~exist(fullPath , 'dir')
    % Folder does not exist so create it.
    mkdir(fullPath );
end

fullPath2 = fullfile(sqpath, 'miniMDS_data');
if ~exist(fullPath2 , 'dir')
    % Folder does not exist so create it.
    mkdir(fullPath2 );
end




fprintf('file name  = %s\n',fullPath,fullPath2);
Resolution = 40000;


 for i =6:13
     fprintf('Processing chromosome %d............ \n', i);
     path1 = strcat(sqpath,'/nij.chr',int2str(i));
     path2 = strcat(fullPath,'/chr',int2str(i),'_list.txt');
     path3 = strcat(fullPath2,'/chr',int2str(i),'_miniMDS.bed');
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
	 
	 
	 
     fprintf('file name  = %s\n', path2);	
     fprintf('file name  = %s\n', path3);		 
	 fileID = fopen(path2,'w');  
	 fileID1 = fopen(path3,'w');  
     for k=1: length(newD) 
          fprintf(fileID,'%d\t%d\t%f',newD(k,:));
          fprintf(fileID,'\n');     
          chr=['chr',int2str(i)];
          fprintf(fileID1,'%6s %12d %12d %6s %12d %12d %9.3f \n',chr,newD(k,1), newD(k,1)+Resolution, chr, newD(k,2), newD(k,2)+Resolution, newD(k,3));  		  
     end 
	 
     fclose(fileID);
	 fclose(fileID1);
 end
   
    
	   
 disp('Completed Successfully');