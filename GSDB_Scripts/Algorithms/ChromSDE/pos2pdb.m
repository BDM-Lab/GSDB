function pos2pdb(posfile)

bigend=10000000000;
Alphabet=char('A'+(1:26)-1)';
[pathstr,fname,ext] = fileparts(posfile);
pdffilename=[pathstr,'/' ,fname, '.pdb'];
fileID = fopen(pdffilename,'w');
posMatrix=load(posfile);
chromBeta=posMatrix(posMatrix(:,1)==-1,2:4);
BACHpos1=posMatrix(posMatrix(:,1)>0,:);
chromlist=chromBeta(:,1);
fprintf(fileID,'AUTHOR    ZHANG ZHIZHUO\n')
startid=0;

s=min(min(BACHpos1(:,2:4)));
r=99/(max(max(BACHpos1(:,2:4)))-min(min(BACHpos1(:,2:4))));
for chrom=chromlist'
     selchrom1=(BACHpos1(:,1)>=chrom*bigend & BACHpos1(:,1)<(chrom+1)*bigend);
     Q3D=BACHpos1(selchrom1,2:4)';
     chromChar=Alphabet(chrom);
for j=1:size(Q3D,2)
       stype2='LIG';
       stype1='C';
       stype3=75.00;
%       fprintf(fileID,'ATOM  %5d %s    %s %c        %8.3f%8.3f%8.3f%6.2f%6.2f    \n',startid+j,stype1,stype2,chromChar,Q3D(1,j),Q3D(2,j),Q3D(3,j),1.00,stype3);    
%        stype2='EDG';
%        stype1='O';
%        stype3=75.00;
       fprintf(fileID,'ATOM  %5d %s    %s %c        %8.3f%8.3f%8.3f%6.2f%6.2f    \n',j,stype1,stype2,chromChar,r*(Q3D(1,j)-s),r*(Q3D(2,j)-s),r*(Q3D(3,j)-s),1.00,stype3);
end
for j=1:(size(Q3D,2)-1)
     fprintf(fileID,'CONECT %4d%4d\n',startid+j,startid+j+1);
end
startid=startid+size(Q3D,2);
end
 fprintf(fileID,'END');

fclose(fileID) ;
   %molviewer(pdffilename)
end