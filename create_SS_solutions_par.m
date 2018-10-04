%this file generates the .mat-file PAR.mat. It contains steady state
%solutions of both models (wild-type and no-allosteric feedback model)
function create_SS_solutions_par
ParSize = 5000; %adjust number of Parameter sets
tspan = 0:1:5000;
par_end = zeros(14,ParSize);
sol_wt = zeros(4,ParSize);
sol_af = zeros(4,ParSize);
opts = odeset('RelTol',1e-06,'AbsTol',1E-6);

parfor i1 = 1:ParSize %parfor is possible
    disp(i1)
    value = 0; %repeat when stability criteria are not met
    while value == 0
        [p,par] = Sample(1,1); %samples one parameter set from case gg = 1
        y0 = par(end-3:end); %initial conditions
        for af = 0:1  %zero refers to the respective mutant
            
            [~,y]  =  ode23s(@(t,c) odemodel(t,c,p,par,af),tspan,y0,opts);
            
            %Parameters;
            beta11  =  par(find(strcmp(p,'beta11')),1);
            beta12  =  par(find(strcmp(p,'beta12')),1);
            beta21  =  par(find(strcmp(p,'beta21')),1);
            beta22  =  par(find(strcmp(p,'beta22')),1);
            K1  =  par(find(strcmp(p,'K1')),1);
            Km  =  par(find(strcmp(p,'Km')),1);
            K2  =  par(find(strcmp(p,'K2')),1);
            alpha1  =  par(find(strcmp(p,'alpha1')),1);
            muemax  =  par(find(strcmp(p,'muemax')),1);
            Kmue  =  par(find(strcmp(p,'Kmue')),1);
            
            %Variables at t(end)
            met1  =  y(end,1);
            met2  =  y(end,2);
            e1  =  y(end,3);
            e2  =  y(end,4);
            
            %Definition of the growth rate
            mue  =  muemax*met2/(met2+Kmue);
            
            %Mass balance
            if af == 1 % wild-type
                dm1dt = beta11*e1*K1/(K1+met2)-beta12*e2*met1/(met1+Km);
            elseif af == 0 %only transcriptional feedback model
                dm1dt = beta11*e1-beta12*e2*met1/(met1+Km);
            end
            dm2dt = beta12*e2*met1/(met1+Km)-alpha1*mue;
            de1dt = beta21*K2/(K2+met2)-e1*mue;
            de2dt = beta22*K2/(K2+met2)-e2*mue;
            
            %check if solutions are in steady state
            F = [dm1dt;dm2dt;de1dt;de2dt];
            if af == 0
                SS_af = sum(abs(F));
            elseif af == 1
                SS_wt = sum(abs(F));
            end
            
            %calculation of the Jacobian Matrix
            if af == 1
                J = [(beta12*e2*met1)/(Km + met1)^2 - (beta12*e2)/(Km + met1),                                                            -(K1*beta11*e1)/(K1 + met2)^2,      (K1*beta11)/(K1 + met2),   -(beta12*met1)/(Km + met1)
                     (beta12*e2)/(Km + met1) - (beta12*e2*met1)/(Km + met1)^2,                     (alpha1*met2*muemax)/(Kmue + met2)^2 - (alpha1*muemax)/(Kmue + met2),                            0,    (beta12*met1)/(Km + met1)
                                                                            0, (e1*met2*muemax)/(Kmue + met2)^2 - (e1*muemax)/(Kmue + met2) - (K2*beta21)/(K2 + met2)^2, -(met2*muemax)/(Kmue + met2),                            0
                                                                            0, (e2*met2*muemax)/(Kmue + met2)^2 - (e2*muemax)/(Kmue + met2) - (K2*beta22)/(K2 + met2)^2,                            0, -(met2*muemax)/(Kmue + met2)];
                
            elseif af == 0
                J = [(beta12*e2*met1)/(Km + met1)^2 - (beta12*e2)/(Km + met1),                                                                                        0,                       beta11,   -(beta12*met1)/(Km + met1)
                     (beta12*e2)/(Km + met1) - (beta12*e2*met1)/(Km + met1)^2,                     (alpha1*met2*muemax)/(Kmue + met2)^2 - (alpha1*muemax)/(Kmue + met2),                            0,    (beta12*met1)/(Km + met1)
                                                                            0, (e1*met2*muemax)/(Kmue + met2)^2 - (e1*muemax)/(Kmue + met2) - (K2*beta21)/(K2 + met2)^2, -(met2*muemax)/(Kmue + met2),                            0
                                                                            0, (e2*met2*muemax)/(Kmue + met2)^2 - (e2*muemax)/(Kmue + met2) - (K2*beta22)/(K2 + met2)^2,                            0, -(met2*muemax)/(Kmue + met2)];
                
            end
            
            %calculate real parts of the Jacobian
            if af == 0
                ew_af = max(real(eig(J)));
                sol_af1 = [met1;met2;e1;e2];
            elseif af == 1
                ew_wt = max(real(eig(J)));
                sol_wt1 = [met1;met2;e1;e2];
            end
            
        end %af
        %check for steady state, eigenvalues (stability) and zero-crossing.
        if SS_wt<1E-08 && SS_af<1E-08 && min(min(y))>0 && ew_wt<-1E-05 && ew_af<-1E-05
            %if true, exit
            value = 1;
        end
    end
    %The steady state solutions are stored in the pre-allocated solution vectors
    sol_af(:,i1) = sol_af1;
    sol_wt(:,i1) = sol_wt1;
end

%steady state solution vectors
%wt
par_end(end-3:end,:) = sol_wt;
par_wt = par_end;
%af
par_end(end-3:end,:) = sol_af;
par_af = par_end;
save('PAR.mat','par_wt','par_af')
end


function dcdt  =  odemodel(~,c,p,par,af,~)

par(end-3:end,1) = c;

%get kinetic parameters
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

%mass balance
if af == 1
    dcdt(1,1) = beta11*e1*K1/(K1+met2)-beta12*e2*met1/(met1+Km);
elseif af == 0
    dcdt(1,1) = beta11*e1-beta12*e2*met1/(met1+Km);
end
dcdt(2,1) = beta12*e2*met1/(met1+Km)-alpha1*mue;
dcdt(3,1) = beta21*K2/(K2+met2)-e1*mue;
dcdt(4,1) = beta22*K2/(K2+met2)-e2*mue;

end

