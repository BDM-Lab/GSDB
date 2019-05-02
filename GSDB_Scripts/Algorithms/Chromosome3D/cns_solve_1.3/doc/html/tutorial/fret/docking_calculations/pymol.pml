load ./sorted/sort_1.pdb, 2
load ./sorted/sort_197.pdb , 2
load ./sorted/sort_215.pdb , 3
load ./sorted/sort_223.pdb , 4


color blue, segid A
color red, segid B
color green, segid C
color lime, segid D

color orange, segid S and resid 140:262
color yellow, segid S and resid 273:418
color white, segid S and resid 262:273

# show sphere, resname DUMY
label resname DUMY, resi

show sticks, ( resid 398 or resid 399 ) and chain S
color red, ( resid 398 or resid 399 ) and chain S

show sticks, ( resid 51 or resid 52 or resid 55 ) and chain C
color lime, ( resid 51 or resid 52 or resid 55 ) and chain C

color red, ( segid S and ( resid 170:177 or resid 231:236 or resid 302:307 or resid 364:369 ))


show cartoon, all

hide lines, all

hide nonbonded, all

# medium FRET 
distance mid, resid 140 and resname DUMY and chain S, resid 20 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 140 and resname DUMY and chain S, resid 139 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 154 and resname DUMY and chain S, resid 20 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 269 and resname DUMY and chain S, resid 20 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 350 and resname DUMY and chain S, resid 139 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 383 and resname DUMY and chain S, resid 20 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 154 and resname DUMY and chain S, resid 197 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 269 and resname DUMY and chain S, resid 197 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 350 and resname DUMY and chain S, resid 197 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 154 and resname DUMY and chain S, resid 139 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 269 and resname DUMY and chain S, resid 139 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 350 and resname DUMY and chain S, resid 20 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 383 and resname DUMY and chain S, resid 139 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 269 and resname DUMY and chain S, resid 76 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 383 and resname DUMY and chain S, resid 197 and resname DUMY and ( chain C or chain D ) 

distance mid, resid 154 and resname DUMY and chain S, resid 61 and resname DUMY and ( chain A ) 

distance mid, resid 154 and resname DUMY and chain S, resid 76 and resname DUMY and ( chain C or chain D )

distance mid, resid 154 and resname DUMY and chain S, resid 396 and resname DUMY and chain S 

distance mid, resid 154 and resname DUMY and chain S, resid 383 and resname DUMY and chain S 

distance mid, resid 189 and resname DUMY and chain S, resid 396 and resname DUMY and chain S 


# low and high FRET

distance low, resid 368 and resname DUMY and chain S, resid 20 and resname DUMY and ( chain C or chain D ) 

distance low, resid 368 and resname DUMY and chain S, resid 139 and resname DUMY and ( chain C or chain D ) 

distance high, resid 383 and resname DUMY and chain S, resid 76 and resname DUMY and ( chain C or chain D ) 

distance low, resid 140 and resname DUMY and chain S, resid 76 and resname DUMY and ( chain C or chain D ) 

distance low, resid 140 and resname DUMY and chain S, resid 197 and resname DUMY and ( chain C or chain D ) 

distance high, resid 350 and resname DUMY and chain S, resid 76 and resname DUMY and ( chain C or chain D ) 

distance low, resid 368 and resname DUMY and chain S, resid 76 and resname DUMY and ( chain C or chain D ) 

distance low, resid 368 and resname DUMY and chain S, resid 197 and resname DUMY and ( chain C or chain D ) 

distance high, resid 269 and resname DUMY and chain S, resid 61 and resname DUMY and ( chain A ) 

distance high, resid 383 and resname DUMY and chain S, resid 61 and resname DUMY and ( chain A ) 

distance low, resid 140 and resname DUMY and chain S, resid 61 and resname DUMY and ( chain A ) 

distance low, resid 368 and resname DUMY and chain S, resid 61 and resname DUMY and ( chain A ) 

distance low, resid 350 and resname DUMY and chain S, resid 61 and resname DUMY and ( chain A ) 

distance low, resid 254 and resname DUMY and chain S, resid 396 and resname DUMY and chain S 


# control distance

distance control, resid 154 and resname DUMY and chain S, resid 174 and resname DUMY and chain S

