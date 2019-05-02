
Path = '/storage/htc/bdm/tosin/GSDB/Data/GSE105697/';
files = dir([Path,'*.h5']);

for i = 1:length(files)
	path='/storage/htc/bdm/tosin/GSDB/Data/GSE105697/'
    name =  files(i).name;
    fprintf('Processing file %d : %s\n', i, name);
    	
	hdfname = strrep(name,'_chromatin_interactions_hg19.h5','');
     
    filename = [Path, name];

    bin_positions = h5read(filename, '/bin_positions');
    bin_positions =  bin_positions';
    chr_bin_range = h5read(filename, '/chr_bin_range');
    chr_bin_range = chr_bin_range';
    chrs = h5read(filename, '/chrs');
    interactions = h5read(filename, '/interactions');
    fprintf(' Resolution = %d\n',bin_positions(1,3) );
    %----------------------------------------------------
    % create a new directory to work
	disp('Done Reading inputs ....');
    path= [path,'Extracted_Data/',hdfname,'/'];
    mkdir(path);
    %----------------------------------------------------
    %======================================================
    % Create a sparse matrix with input
    %======================================================
    n = length(interactions);


    genomename = [hdfname,'_All_Genome.txt'];
    name = [path,genomename];
    fid = fopen(name,'w');
    for i = 1:length(interactions)
        for j = i:length(interactions)
            Sparse = [ i, j, interactions(i,j)];
            fprintf(fid,'%d \t %d\t %.1f\n',Sparse);
        end     
    end

    fclose(fid);


    % write sequence length to file
    sequencelen =[];

    %======================================================
    % Create each Chromosome file
    %======================================================
    for i = 1:length(chr_bin_range)    
       chromosomename = ['chr',int2str(i),'.txt'];   
        if (i > 22)
            if (i==23)
               chromosomename = 'chrX.txt';
            elseif(i == 24)
                chromosomename = 'chrY.txt';
            else
                chromosomename = 'chrM.txt';
            end

        end 
         name = [path, chromosomename];
         fid = fopen(name,'w');

        range = chr_bin_range(i,:);
        Start = range(1,1) + 1; End = range(1,2) + 1;

        sequencelen = [sequencelen, End - Start + 1];

        for j = Start:End
            for k =  Start:End
                 row_data = [bin_positions(j,2),bin_positions(k,2), interactions(j,k)]; % start start interaction
                 fprintf(fid,'%d \t %d\t %.1f\n',row_data);
            end
        end
        fclose(fid);
    end

    %======================================================
    % Create sequence length
    %======================================================
     lengthname = [hdfname,'_chrom_sequence_length.txt'];
     name = [path,  lengthname];
     dlmwrite(name, sequencelen); % Get chromosome length

    %======================================================
    % Create mapping matrix 
    %======================================================
    mappingname = [hdfname,'_mapping','.txt'];
     name = [path, mappingname];
     fid = fopen(name,'w');
    for i = 1: length(bin_positions)
        order = [i, bin_positions(i,1)+1, bin_positions(i,2),bin_positions(i,3)];
        fprintf(fid,'%d \t %d \t %d \t %d\n', order);
    end
     fclose(fid);

     disp('Completed Successfully!!!!');
end