file = int2str(22);
filename= strcat(file,'.mat');
c= load(filename);
ds = struct2cell(c); % Convert Struct to Mat
contact= cell2mat(ds); 
normalization_3DEM; %normalzation
name = 'best_structure_chr22.pdb';
XYZ = extract_coordinates(name); %Get coordinates
IF_threshold =0.65;
[ Contact_Satisfaction_Score, NonContact_Satisfaction_Score ...
    , IF_Satisfaction_Score]= Scoring_function_3DEM(contact,XYZ,IF_threshold)