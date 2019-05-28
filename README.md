
## GSDB : Genome Structure Database

----------

#### Bioinformatics, Data Mining, Machine Learning (BDM) Laboratory, 
#### University of Missouri, Columbia MO 65211

----------

#### Contact: <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Jianlin Cheng, PhD <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Department of Computer Science <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; University of Missouri, Columbia <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Email: chengji@missouri.edu


GSDB is a database of Hi-C data chromosome and genome structures. Our goal is that this database will enable the exploration of the dynamic architecture of the different Hi-C 3D structure in a variety of cells and tissues.
<br/><br/>
Over 50,000 structures from 12 start-of-the-art Hi-C data structure prediction algorithms for 32 Hi-C datasets each containing varying resolutions.



## 1. Content of folders:
- **src**: source code for the website
- **GSDB_Scripts**: Contains the Algorithms and the scripts used for data extraction, data normalization, and 3D structure generation.
- **Database Data Info.xlsx**: Contains more information about the Hi-C data and the 3D structure prediction tools
- **ID_Generator.jar**: Java executable file to create the GSDB ID. To execute type in terminal: java -jar ID_Generator.jar

## 2. Algorithms Input ##

For each algorithm, we have described the contact matrix input format they accept and the input file name extension/suffix used for the 3D structure Construction

| Algorithm|  Input Format	| GSDB Input filename suffix |
| --- | --- |--- |
|LorDG | 3-column Matrix | _list.txt |
|3DMax| 3-column Matrix| _list.txt|
|MOGEN| 3-column Matrix| _list.txt|
|Pastis|3-column Matrix(bin1,bin2,IF), and mapping coordinate(chr, start_pos,end_pos, bin) | .n_contact, .cbins|
|Chromosome3D| n x n Square Matrix| _matrix.txt|
|HSA| 2-column bin positon with n x n Square Matrix|_HSA.txt|
|miniMDS|Chromosome,positon,IF(chr,start_pos1,end_pos1,chr, start_pos2,end_pos2,IF)|.bed|
|ShRec3D|n x n Square Matrix| _matrix.txt|
|GEM|2-column bin positon with n x n Square Matrix|_HSA.txt|
|ChromSDE|3-column Matrix, and mapping coordinate|.n_contact, .cbins|
|SIMBA3D	|n x n Square Matrix|.npy|
|InfMod3DGen|n x n Square Matrix| _matrix.txt|

## 3. Disclaimer ##

The executable software and the source code of is distributed free of charge as it is to any non-commercial users. The authors hold no liabilities to 
the performance of the program.


