! To compute the IR frequency of the structure
%mem=24GB
%oldchk=optim_fin.chk
%chk=freq.chk
#p UB3LYP/6-31+G(2df,p) opt=(CalcAll,maxcycle=200,maxstep=20) scrf=(pcm,solvent=THF) pop=reg freq=noraman

Structure of plone0

-3  1
P	-0.000260   -0.000780    0.609798
S        0.029369    1.995065   -0.190366
S       -1.745495   -0.971852   -0.190714
S        1.716370   -1.022482   -0.190606
 
