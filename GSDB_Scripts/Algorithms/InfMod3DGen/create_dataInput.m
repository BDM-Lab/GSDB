list_data=dlmread(inputdata);
list_data(:,1)=list_data(:,1)+1;
list_data(:,2)=list_data(:,2)+1;
N=length(list_data);
column1=repmat(ChrNo,N,1);
column2=list_data(:,1);
column3=repmat(ChrNo,N,1);
column4=list_data(:,2);
column5=list_data(:,3);
newD=[column1 column2 column3 column4 column5];
fileID = fopen(filename,'w');  
 for k=1: N 
      fprintf(fileID,'%d\t%d\t%d\t%d\t%f',newD(k,:));
      fprintf(fileID,'\n');  
 end 
 fclose(fileID);

%Use NewD as Input Data

inputdata=filename;