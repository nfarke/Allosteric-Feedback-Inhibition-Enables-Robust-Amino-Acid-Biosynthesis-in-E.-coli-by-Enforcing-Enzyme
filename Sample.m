function [p,par]=Sample(ParSize,gg)

p={'beta11' 'beta12' 'beta21' 'beta22'   'K1' 'Km' 'K2'  'alpha1' 'muemax' 'Kmue'   'met1' 'met2' 'e1' 'e2'};

%parameters average flux
if gg == 1 %called by create_SS_solutions_par -> for figure4b
    par_lb=[930     930	    0.0085/100*1	0.0085/100 0.01	 0.01	0.01	 	86	0.01  1E-5	 0.01	0.01	1.00E-05	1.00E-05]';
    par_ub=[4140    4140   0.0085/100*10  0.0085/10    1      1      1  	86	0.01  1E-5   0.01	0.01    1.00E-05	1.00E-05]';
elseif gg == 0 %called by create_SS_solutions_par2 -> for figure4c
    par_lb=[930     930	    0.0085/100*1	0.0085/100  0.0001	 0.01	0.0001	  	86	0.01  1E-5	 0.01	0.01	1.00E-05	1.00E-05]';
    par_ub=[4140    4140    0.0085/100*10  0.0085/10   10        1      10      86	0.01  1E-5   0.01	0.01    1.00E-05	1.00E-05]';
end

par=zeros(length(par_lb),ParSize);
%random sampling
for i=1:ParSize
    par(:,i)= 10.^((log10(par_lb)-log10(par_ub)).*rand(length(par_lb), 1)+log10(par_ub));
end

end