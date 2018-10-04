function create_SS_solutions_par2
%this function generations steady state solutions and parameter sets for
%the continuation method (PAR.mat)

ParSize = 5000;
tspan = 0:1:5000;
%pre-allocate solution vectors
sol_wt = zeros(4,ParSize);
ew = zeros(1,ParSize);
gg = 0;
[p,~] = Sample(ParSize,gg); %parameter sampling
PAR = zeros(14,ParSize);

opts = odeset('RelTol',1e-07,'AbsTol',1E-7);

parfor i1 = 1:ParSize
    value = 0;        
    while value == 0
        [~,par] = Sample(1,gg);
        y0 = par(end-3:end);
        disp(i1);
        
        [t,y]  =  ode23s(@(t,c) odemodel(t,c,p,par,i1),tspan,y0,opts);
        
        beta11 = par(find(strcmp(p,'beta11')),1);
        beta12 = par(find(strcmp(p,'beta12')),1);
        beta21 = par(find(strcmp(p,'beta21')),1);
        beta22 = par(find(strcmp(p,'beta22')),1);
        K1 = par(find(strcmp(p,'K1')),1);
        Km = par(find(strcmp(p,'Km')),1);
        K2 = par(find(strcmp(p,'K2')),1);
        alpha1 = par(find(strcmp(p,'alpha1')),1);
        muemax = par(find(strcmp(p,'muemax')),1);
        Kmue = par(find(strcmp(p,'Kmue')),1);
        
        met1 = y(end,1);
        met2 = y(end,2);
        e1 = y(end,3);
        e2 = y(end,4);
        
        mue = muemax*met2/(met2+Kmue);
        
        %check, if solutions are in steady state at t_end. 4 mass bilance
        %equations.
        dmet1dt = beta11*e1*K1/(K1+met2)-beta12*e2*met1/(met1+Km);
        dmet2dt = beta12*e2*met1/(met1+Km)-alpha1*mue;
        de1dt = beta21*K2/(K2+met2)-e1*mue;
        de2dt = beta22*K2/(K2+met2)-e2*mue;
        
        F = [dmet1dt;dmet2dt;de1dt;de2dt];
        SS = sum(abs(F));
        
        %Calculate Jacobian Matrix J for each condition and parameterisation 
        %to make sure the steady state is sufficiently stable.
        J = [ (beta12*e2*met1)/(Km + met1)^2 - (beta12*e2)/(Km + met1),                                                            -(K1*beta11*e1)/(K1 + met2)^2,      (K1*beta11)/(K1 + met2),   -(beta12*met1)/(Km + met1)
              (beta12*e2)/(Km + met1) - (beta12*e2*met1)/(Km + met1)^2,                     (alpha1*met2*muemax)/(Kmue + met2)^2 - (alpha1*muemax)/(Kmue + met2),                            0,    (beta12*met1)/(Km + met1)
                                                                     0, (e1*met2*muemax)/(Kmue + met2)^2 - (e1*muemax)/(Kmue + met2) - (K2*beta21)/(K2 + met2)^2, -(met2*muemax)/(Kmue + met2),                            0
                                                                     0, (e2*met2*muemax)/(Kmue + met2)^2 - (e2*muemax)/(Kmue + met2) - (K2*beta22)/(K2 + met2)^2,                            0, -(met2*muemax)/(Kmue + met2)];
        
        %calculate eigenvalues
        ew(i1) = max(real(eig(J)));
        
        %check for steady state, eigenvalues (stability), code duration (needs
        %to run to the end) and zero-crossing. If satisfied the steady state
        %solutions are stored.
        if SS<1E-08 && ew(i1)<-1E-5 && t(end) == tspan(end)  && min(min(y))>0
            PAR(:,i1) = par;
            sol_wt(:,i1) = [met1;met2;e1;e2];
            value = 1;
        end
        
    end
    
end

PAR(end-3:end,:) = sol_wt;
par = PAR(:,:);
%ratio of feedback parameters
ratio = par(find(strcmp(p,'K2')),:)./par(find(strcmp(p,'K1')),:);
save('PAR2.mat','par','ratio');
end



function dcdt  =  odemodel(t,c,p,par,i1,opts)

par(end-3:end,1) = c;
beta11 = par(find(strcmp(p,'beta11')),1);
beta12 = par(find(strcmp(p,'beta12')),1);
beta21 = par(find(strcmp(p,'beta21')),1);
beta22 = par(find(strcmp(p,'beta22')),1);
K1 = par(find(strcmp(p,'K1')),1);
Km = par(find(strcmp(p,'Km')),1);
K2 = par(find(strcmp(p,'K2')),1);
alpha1 = par(find(strcmp(p,'alpha1')),1);
muemax = par(find(strcmp(p,'muemax')),1);
Kmue = par(find(strcmp(p,'Kmue')),1);
e1 = par(find(strcmp(p,'e1')),1);
e2 = par(find(strcmp(p,'e2')),1);
met2 = par(find(strcmp(p,'met2')),1);
met1 = par(find(strcmp(p,'met1')),1);
mue = muemax*met2/(met2+Kmue);

dcdt(1,1) = beta11*e1*K1/(K1+met2)-beta12*e2*met1/(met1+Km);
dcdt(2,1) = beta12*e2*met1/(met1+Km)-alpha1*mue;
dcdt(3,1) = beta21*K2/(K2+met2)-e1*mue;
dcdt(4,1) = beta22*K2/(K2+met2)-e2*mue;

end

