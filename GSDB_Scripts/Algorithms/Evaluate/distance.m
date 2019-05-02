% Extract the non NAN 
ntfd_index=[];
Y=[];
for i=1:length(XYZ)
	if (isnan(XYZ(i,1)))
		ntfd_index=[ntfd_index, i];
		continue;
	else
		 Y= [Y;XYZ(i,:)];    
	end       
end

n=length(Y);
XYZ=zeros(n,3);
XYZ(:,1)=Y(:,1);
XYZ(:,2)=Y(:,2);
XYZ(:,3)=Y(:,3);