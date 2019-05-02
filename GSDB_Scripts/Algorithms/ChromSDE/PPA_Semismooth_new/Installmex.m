%%******************************************************
%% Run this script in Matlab command window
%%
%%******************************************************

function Installmex

computer_model = computer;
matlabversion = sscanf(version,'%f');
matlabversion = matlabversion(1);
%%
if strcmp(computer_model,'PCWIN')
    %       str1 = ['''',matlabroot,'\extern\lib\win32\lcc\libmwlapack.lib''  '];
    %       str2 = ['''',matlabroot,'\extern\lib\win32\lcc\libmwblas.lib''  '];
    str1 = ['''',matlabroot,'\extern\lib\win32\microsoft\libmwlapack.lib''  '];
    str2 = ['''',matlabroot,'\extern\lib\win32\microsoft\libmwblas.lib''  '];
    libstr = [str1,str2];
elseif strcmp(computer_model,'PCWIN64')
    %str1 = ['''',matlabroot,'\extern\lib\win64\lcc\libmwlapack.lib''  '];
    %str2 = ['''',matlabroot,'\extern\lib\win64\lcc\libmwblas.lib''  '];
    str1 = ['''',matlabroot,'\extern\lib\win64\microsoft\libmwlapack.lib''  '];
    str2 = ['''',matlabroot,'\extern\lib\win64\microsoft\libmwblas.lib''  '];
    libstr = [str1,str2];
else
    libstr = '-lmwlapack -lmwblas  ';
end
mexcmd = 'mex -O  -largeArrayDims  -output ';
%%

fsp = filesep;

curdir = pwd;
fprintf(' current directory is:  %s\n',curdir);
%%
%% generate mex files in mexfun
%%
clear fname

src = [curdir,fsp,'solver/mexfun'];
eval(['cd ','solver/mexfun']);
fprintf ('\n Now compiling the mexFunctions in:\n');
fprintf (' %s\n',src);
%%
if (ispc)
    cmd([mexcmd, 'mexeig mexeigwin.c','  ',libstr]);
    cmd([mexcmd, 'mexsvd mexsvdwin.c','  ',libstr]);
    cmd([mexcmd, 'mexAV2 mexAV2win.c','  ',libstr]);
else
    cmd([mexcmd, 'mexeig mexeig.c','  ',libstr]);
    cmd([mexcmd, 'mexsvd mexsvd.c','  ',libstr]);
    cmd([mexcmd, 'mexAV2 mexAV2.c','  ',libstr]);
end

%%
fname{1} = 'mexbwsolve';
fname{2} = 'mexfwsolve';
fname{3} = 'mexmat';
fname{4} = 'mexsmat';
fname{5} = 'mexsvec';
fname{6} = 'mextriang';
fname{7} = 'mextriangsp';
fname{8} = 'mexScalarxMatrix';
fname{9} = 'mexhouse';
fname{10} = 'mexMatvec';
fname{11} = 'mexProjOmega';
fname{12} = 'mexspconvert';

for k = 1: length(fname)
    cmd([mexcmd,'  ',fname{k},'  ',fname{k},'.c','  ',libstr]);
end
fprintf ('\n Compilation of mexFunctions completed.\n');
cd ..
cd ..
%%***********************************************
function cmd(s)

fprintf(' %s\n',s);
eval(s);
%%***********************************************
