   %Extract the cbins from matrix
	Data1=contact(:,1:2);
	contact=contact(:,3:end);
	len = length(contact);
	
	%create	c contact
    count = 0;
    ncontact = [];
    for i = 1:len    
        for j = i:len
            count = count + 1;
            c=[i,j,contact(i,j)];
            ncontact = [ncontact;c ];
    %       
        end
    end

	
	%create c contact
	% Note: Reordered the input cbins so  that no reordering required in the construction code
	chrome=repmat(nn,len,1);
	Len = [1:len]'	;
	Data = [chrome Data1 Len];
	
	% Save the Data to specified path
	if (nn==24)	
	 nn='Y';
	end
	if (nn==25)	
	 nn='M';
	end
	
    pathname =strcat(output_path,'/','chr',num2str(nn),'.cbins');
    fileID = fopen(pathname,'w');    
	 for j = 1:length(Data(:,1))     
		fprintf(fileID, '%d \t %d \t %d \t %d ',Data(j,1:end));  
	  fprintf(fileID, '\n');   
	 end          
     fclose(fileID);
	
    fragname = strcat(output_path,'/','chr',num2str(nn),'.n_contact');
    dlmwrite(fragname,ncontact,'delimiter','\t');