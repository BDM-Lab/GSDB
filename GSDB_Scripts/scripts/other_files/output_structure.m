
XYZ=XYZ*10;
%output pdb file
data.X=XYZ(:,1);
data.Y=XYZ(:,2);
data.Z=XYZ(:,3);
[pathstr,name,ext] = fileparts(filename);

data.outfile =  ['Structure_5kb/',name,'.pdb']; %output directory.
if exist(data.outfile, 'file')==2   %delet file if exists.
  delete(data.outfile);
end

mat2pdb(data); % Converts the mat XYZ coordinates to PDB format.
