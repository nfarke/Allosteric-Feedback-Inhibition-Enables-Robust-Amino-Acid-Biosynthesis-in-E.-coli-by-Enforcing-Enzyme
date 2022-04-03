Allosteric Feedback Inhibition Enables Robust Amino Acid Biosynthesis in E. coli by Enforcing Enzyme Overabundance
DOI: 10.1016/j.cels.2018.12.005

I built an ODE-based dynamic mathematical model of the amino acid biosynthesis in E. coli. The model is a simplified
representation of the amino acid pathways with two enzymatic conversions to the amino acid. The amino acid allosterically
inhibits the first reaction and both enzyme production rates (enzyme level regulation). Rate laws and regulation were
modelled according to michaelis-menten type kinetics.

We then created three different models, the wild-type model (WT), the no-allosteric feedback model were the allosteric
feedback was removed and the no-enyzme-level regulation model without enzyme level regulation.
We then sampled 5000 parameter sets from a logarithmic distribution and all models shared the same steady state flux.
All models satisfied a linear stability constraint.
The different mutants showed different enzyme and metabolite levels which matched the experimental observations.

We then used a parameter continuation method to simulate a bottleneck in the pathways for 5000 parameter sets in 
the WT model. The result was a trade-off between robustness and efficiency(enzyme amount)




-------------------------------------------------------------------------------------------------------------

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
