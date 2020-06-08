function Figure4c
load PAR2.mat
Pert = 500; %perturbation degree 500 fold in this case
NoSteps = Pert*10; %number of steps
varpert = {'beta22'}; %define perturbation variable, here it is beta22
varpert1 = sym(varpert);


p = {'beta11' 'beta12' 'beta21' 'beta22'   'K1' 'Km' 'K2'  'alpha1' 'muemax' 'Kmue'   'met1' 'met2' 'e1' 'e2'};
syms 'beta11' 'beta12' 'beta21' 'beta22' 'K1' 'Km' 'K2' 'alpha1' 'muemax' 'Kmue' 'met1' 'met2' 'e1' 'e2'
EnsembleSize=size(par,2);
par1 = par;
psym = sym(p);

%set up rate equations
r1 = beta11*e1*K1/(K1+met2);
r2 = beta12*e2*met1/(met1+Km);
mue = muemax*met2/(met2+Kmue);
r3 = alpha1*mue;
rTR1 = beta21*(K2/(K2+met2));
rTR2 = beta22*(K2/(K2+met2));
dil1 = e1*mue;
dil2 = e2*mue;

%define rate vector r, stoichiometric matrix S and mass balance F = S*r
r = [r1;r2;r3;rTR1;rTR2;dil1;dil2];
S = [1 -1 0 0 0 0 0; 0 1 -1 0 0 0 0; 0 0 0 1 0 -1 0; 0 0 0 0 1 0 -1];
F = S*r; 

%set up derivatives as matlab function for enhanced speed
DFDX = matlabFunction(jacobian(F,[met1,met2,e1,e2]),'vars',{psym});
DFDP = matlabFunction(jacobian(F,varpert1),'vars',{psym});

%pre-allocate Result matrices
ModelResults = repmat({NaN(1,NoSteps+1)},[EnsembleSize 1]);
conc=repmat({NaN(1001,5)},[EnsembleSize 1]);

parfor n = 1:EnsembleSize
    disp(n)    
    par = par1(:,n);
    Varini = par(find(strcmp(p,varpert))); %initial state of pert variable
    Varmax = par(find(strcmp(p,varpert)))/Pert; %state at maximum level of perturbation
    InitialState = [par(end-3) par(end-2) par(end-1) par(end)];
    
    options = odeset('Events',@event_function,'RelTol', 1e-7, 'AbsTol', 1e-7);
    TimeIn=clock;
    [t, conc1] = ode23s(@dxFunc,0:1/NoSteps:1,InitialState,options,psym,par,Varini,Varmax,DFDX,DFDP,TimeIn,varpert);
    Var = repmat(Varini,1,length(t)) + repmat(t',length(Varini),1).*repmat(Varmax-Varini,1,length(t));
    Var = Var./Varini; %normalize results to initial state
    ModelResults(n,:) = {Var}; %distance
    conc(n,:) = {conc1}; %vector containing variables
end

for iii=1:EnsembleSize
    pert_final(iii,1) = 1-min(ModelResults{iii,1});
end

cx = exp_colormap('blue-orange',EnsembleSize/100); %calls file exp_colormap.m
[ratio1,I] = sort(ratio,'ascend');
e_start = sum(par1(end-1:end,I));
pert_final1 = pert_final(I);

ratio1_reshaped=reshape(ratio(1:length(ratio1)),[100 EnsembleSize/100]);
e_start_reshaped=reshape(e_start(1:length(e_start)),[100 EnsembleSize/100]);
pert_final_reshaped=reshape(pert_final1(1:length(pert_final1)),[100 EnsembleSize/100]);

mte=0.17;% maximum theoretical enzyme (mte) amount (e1+e2, taking the
%biggest beta values and assuming complete derepression

%create figure
for iv = 1:size(ratio1_reshaped,2)
    scatter(e_start_reshaped(:,iv)/mte*100,pert_final_reshaped(:,iv)*100,10,cx(iv,:),'filled');
    hold on
end
index=1:EnsembleSize/10.005:EnsembleSize;
set(gca,'xdir','reverse');
colormap(cx);
set(gca,'YTick',[0 25 50 75 100],'YTickLabels',[0 25 50 75 100])
colorbar('yticklabel',round(log10(ratio1(round(index))),1));
set(gca,'xscale','log');
xlim([round(min(e_start)/mte) 100])
ylabel('Robustness (%)');
xlabel('Enzyme level (%)');
end



function dx= dxFunc(t,x,psym,par,Varini,Varfinal,DFDX,DFDP,~,varpert)
par(end-3:end) = x;
U = Varini + t.*(Varfinal-Varini);
par(find(ismember(psym,varpert))) = U;
XJac = DFDX(par');
DFDP1 = DFDP(par');
dx = -XJac\DFDP1*(Varfinal-Varini);
end

function [value,isterminal,direction] = event_function(t,x,psym,par,Varini,Varfinal,DFDX,~,TimeIn,varpert)
par(end-3:end)=x;
U = Varini + t.*(Varfinal-Varini);
par(find(ismember(psym,varpert))) = U;
XJac = DFDX(par');
value = max(real(eig(XJac))); % when value = 0, an event is triggered
if max(real(eig(XJac)))> -1E-05 || any(par<0) %|| etime(clock,TimeIn)> 1  
    value = 0;
end
isterminal = 1; % terminate after the first event
direction = 0;  % get all the zeros
end
