
% score = [];
for Chr =2:2
    
    Res=5000;

    fprintf('Chromosome %d....\n',Chr);
	Output_path='/storage/htc/bdm/tosin/High_Res_3DGenome/Outputs/miniMDS/5kb/';
    filename=[Output_path,'chr',num2str(Chr),'_5kb_miniMDS_structure.tsv'];  
    M = dlmread(filename, '', 3);

    disp('Reading observed data.......')
    path=['/storage/htc/bdm/tosin/High_Res_3DGenome/GM12878_5kb_10kb/'];
    fname=[path,'chr',num2str(Chr),'_5kb.KRnormFULL'];
    Data=dlmread(fname);

    start=min(Data(:,1));

     D = unique(Data(:,1:2));
     D = sort(D);
     res=repmat(Res,length(D),1);
     newD = [D,D+res];
     
    % Deleting data
    disp('Deleting data.......')
    clear Data;
    
    disp('Creating new coordinate data.......')
    NewD=[];
    for k=1:length(M)
        End = start + Res;
        range = [];
        NewD = [NewD; start End M(k,2:end)];
        start=End;
    end

    % Extract the non NAN 
    Y=[];
    for i=1:length(NewD)
        if (isnan(NewD(i,3)))
            continue;
        else
             Y= [Y;NewD(i,:)];    
        end       
    end


    disp('Loading SQ data.......')
    SQ_path = '/storage/htc/bdm/tosin/High_Res_3DGenome/GM12878_5kb_10kb_SquareMatrix/';
    filenm=[SQ_path,'chr',num2str(Chr),'_5kb.KRnormSQUARE.txt'];
    %load HSA data
    SQ = dlmread(filenm);

    %Get the position of the last in Square MAtrix
    pos1=NewD(end,1);
    ind = find(newD(:,1)==pos1);

    % Keep old Sqaure matrix
    old_SQ = SQ;

    % cut the matrix at the ind index
    SQ=SQ(1:ind,1:ind);

    %find the matrix NewD on the Square Matrix Data
    fd_index=[];
    ntfd_index=[];
    for m=1:ind
        cur=newD(m,1);
        fd=find(Y(:,1)==cur);
        if (~isempty(fd))
            fd_index=[fd_index; fd]; 
        else
            ntfd_index=[ntfd_index m];
        end
    end


    disp('Removing unused data.............');
    clear old_SQ newD NewD M  Data
    %Remove the rows not on the cordinate emove a row from a matrix
    disp('Now processing SQ data removing missing data.............');
       
    SQ(ntfd_index,:)=[];

    SQ(:,ntfd_index)=[];

    % the matrix and the xyz cordinate should be of same length

    n=length(Y);
    XYZ=zeros(n,3);
    XYZ(:,1)=Y(:,3);
    XYZ(:,2)=Y(:,4);
    XYZ(:,3)=Y(:,5);

    % Check if non_NAN and square matrix are same
    if (n~=length(SQ))
        disp('The Matrices are not of the same length');
    end
  
     disp('Computing sturcture distance.............');
     [dist] = StructureDistance(XYZ,n); % distance from structure
    % 
    %output pdb and image  
    output_structure; 
    %convert IF to distancerows=length(HSA(:,1));
    rows=length(SQ(:,1));
    col=length(SQ(1,:));
    Data=zeros(rows,col);

    % converting contact frequency to distance 
    for j = 1:n
        for k = 1:n
              if(SQ(j,k)~=0)
                 Data(j,k)=1./(SQ(j,k).^4);
              end
        end
    end

    clear SQ 
    [RHO,RHO1] = Spearman_corr(Data,dist);  
    score=[score;RHO];
    fprintf('The spearman correlation for chromsome %d = %f\n',Chr, RHO );
end


