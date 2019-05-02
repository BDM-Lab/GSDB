%format data position properly

path = '/storage/htc/bdm/tosin/High_Res_3DGenome/1kb_normalized/';


 for i = 1:23
   
     fprintf('Processing chromosome %d............ \n', i);
     fpath = strcat(path,'chr',int2str(i));
     datapath = [fpath,'_1kb.KRnormFULL'];
     Data = dlmread(datapath);
     
    bpath='/storage/htc/bdm/tosin/High_Res_3DGenome/1kb_normalized_formatted/';
	opath = strcat(bpath,'chr',int2str(i));
    newpath = [opath,'_1kb.KRnormFULL'];
	 fileID = fopen(newpath,'w');
	 for k=1:length(Data)      
		 fprintf(fileID,'%8d %9d %9.3f\n', Data(k,1), Data(k,2), Data(k,3));        
	 end 
	 fclose(fileID);
 end