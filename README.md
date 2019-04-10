
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
- **GSDB_Scripts**: All the scripts and parameter settings used for 3D structure generation for each algortihm
- **Algorithms**:  The algorithms and parameter settings used for 3D sructure prediction


## 2. Algorithms Input ##

For each algorithm, we have described the contact matrix input format they accept and the input file name extension/suffix used for the 3D structure Construction

| Algorithm|  Input Format	| GSDB Input filename suffix |
| --- | --- |--- |
|LorDG | 3-column Matrix | |
|3DMax| 3-column Matrix||
|MOGEN| 3-column Matrix||
|Pastis|||
|Chromosome3D| n x n Square Matrix||
|HSA| 2-column Bin posiiton with n x n Square Matrix||
|miniMDS|||
|ShRec3D|||
|GEM|||
|ChromSDE|||
|SIMBA3D	|||
|InfMod3DGen|||

## 3. Disclaimer ##

The executable software and the source code of is distributed free of charge as it is to any non-commercial users. The authors hold no liabilities to 
the performance of the program.


