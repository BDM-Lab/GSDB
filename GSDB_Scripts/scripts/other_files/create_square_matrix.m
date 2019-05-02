%% Create a matrix data for the sparse dataset
path = '/storage/htc/bdm/tosin/High_Res_3DGenome/';

Length = [];

 for i = 17:19
     fprintf('Processing chromosome %d............ \n', i);
     fpath = strcat(path,'1kb_normalized/chr',int2str(i));
     datapath = [fpath,'_1kb.KRnormFULL'];
	 fprintf('Processing chromosome file %s............ \n', datapath);
	 
     Data = dlmread(datapath);     
     D = unique(Data(:,1:2));
     ind=sort(D);
     n=length(ind);
     Mat = zeros(n);
     
     for k=1:length(Data)
         i1=find(ind==Data(k,1));
         i2=find(ind==Data(k,2));
         Mat(i1,i2)=Data(k,3);
         Mat(i2,i1)=Data(k,3);
     end
     
     fpath = strcat(path,'1kb_Square/','chr',int2str(i),'_1kb.KRnormSQUARE.txt');
     dlmwrite(fpath,Mat,' ');
 end