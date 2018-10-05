This file contains instructions to replicate the figures from the paper. The code was written in MatlabR2017b

Figure4b:
requires PAR.mat, Sample.m (0) and Figure4b.m
- running Figure4b.m reproduces the raw boxplot from Figure4b

Figure4c
requires PAR2.mat, Sample.m (1), exp_colormap.m and Figure4c.m and the symbolic toolbox
- Figure4c.m can be run with or without parfor. If the parallel toolbox is not available change parfor to for


The .mat-file PAR is generated from create_SS_solutions_par.m
The .mat-file PAR2 is generated from create_SS_solutions_par2.m
Both are compatible with parfor (parallel computing toolbox)
