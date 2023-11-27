clc
close all
clear all

OCN = build_OCN("OCN_A.mat",30*10000*10000);

%%
par = common_parameters(); 





rng('default')
setup = build_setup(OCN,par,33*800*1000,'seed',3108);
rng('default')

setup.T=diag(min(setup.par.U*setup.H,setup.chi.*setup.F));

eqi_RA = zeros(OCN.nNodes,4);
for nn = 1:OCN.nNodes
    tmp = max(find_EE(setup.par,setup.H(nn),setup.S(nn),setup.A(nn),...
        setup.chi(nn),setup.par.epsilon(nn),setup.par.theta(nn),...
        setup.par.xi(nn),setup.T(nn,nn)));
    if  isempty(tmp) || tmp(1)==0
        eqi_RA(nn,:)=NaN;
    else
        eqi_RA(nn,:)=tmp;
    end
end

eqi_vr = int8(~isnan(eqi_RA(:,1))); eqi_vr(eqi_vr==0)= NaN;

WH = sum(eqi_RA(:,1).*setup.H,'omitmissing')/sum(setup.H.*double(eqi_vr),'omitmissing')
YS = sum(eqi_RA(:,3).*setup.S,'omitmissing')/sum(setup.S.*double(eqi_vr),'omitmissing')
IF = sum(eqi_RA(:,4).*setup.F,'omitmissing')/sum(setup.F.*double(eqi_vr),'omitmissing')

sum(~isnan(eqi_RA(:,1)))

sum(setup.sigma>0)


%%

%%

N_iter=1000; % number of iterations for the simulated annealing algorithm
cooling_rate=1/1200; 

dF_MC = zeros(N_iter,1); c_MC = zeros(N_iter,1);
thetaC_MC = zeros(N_iter,1); betaE_MC = zeros(N_iter,1);
EVAL = zeros(N_iter,1);
T_SA_acc = zeros(N_iter,1);
if eqi(1) == 0
    val = -99;
else
    val = -sum([log10(eqi(1))-4 log10(eqi(3))-log10(0.003) log10(eqi(4))-log10(0.269*3)].^2);
end

c_MC(1)=setup.par.c; 
thetaC_MC(1)=setup.par.theta_C;  betaE_MC(1)=setup.par.beta_E;

EVAL(1) = val;

val_current = val;
Temp_SA=@(x) exp(-x*cooling_rate);
T_SA_acc(1) = Temp_SA(0);

indK=0; OldTime=0; OldK=0;
tic
for indI = 1:N_iter

    setup.par.c    = TruncNormRnd(setup.par.c,5e-8,1e-8,1e-6);
    setup.par.U    = setup.par.c*setup.F(node);
    setup.par.theta_C = TruncNormRnd(setup.par.theta_C, 5e-9, 1e-8, 1e-6);
    setup.par.beta_E = TruncNormRnd(setup.par.beta_E, 5e-9, 1e-9, 1e-7);


    setup.KF = par.dF*OCN.SC_RiverLength.*OCN.SC_RiverWidth;   
    setup.chi = setup.par.c*setup.H;
    setup.F = find_fish_equilibrium(setup.KF,setup.W,setup.par,setup.chi);
    
    eqi = max(find_EE(setup.par, setup.chi(node),...
        setup.H(node),setup.S(node),setup.F(node),setup.A(node),...
        setup.par.theta_C,setup.par.xi(node),setup.T(node,node)));
    
    if eqi(1) == 0
        val = -99;
    else
        val = -sum([log10(eqi(1))-4 log10(eqi(3))-log10(0.003) log10(eqi(4))-log10(0.269*3)].^2);
    end
    
    c_MC(1)=setup.par.c; 
    thetaC_MC(1)=setup.par.theta_C;  betaE_MC(1)=setup.par.beta_E;
    
    EVAL(1) = val;
    % accept or reject new parameter set
    if  val>val_current || rand < exp((val-val_current)/Temp_SA(indI))
            indK=indK+1;
            c_MC(indK)=setup.par.c; 
            thetaC_MC(indK)=setup.par.theta_C;  betaE_MC(indK)=setup.par.beta_E;
            EVAL(indK) = val;
            T_SA_acc(indK) = Temp_SA(indI);
    end
    
    % print information on MC
    if mod(indI,10)==0
        NewTime=toc; 
        fprintf('\n')
        fprintf('%0.2f %% COMPLETE\n',indI/N_iter*100)       
        fprintf('Iterations: %d    Chain length: %d   Acceptance rate %0.2f%%\n',indI,indK,indK/indI*100)
        fprintf('Total elapsed time = %.1f s    Average speed = %0.1f iter/s\n',NewTime,indI/NewTime)
        fprintf('Local accept. rate = %.2f%%    Local speed = %0.1f iter/s\n',(indK-OldK)/1000*100,1000/(NewTime-OldTime))
        fprintf('EVAL = %0.4f    Temperature = %0.4f\n',EVAL(indK),Temp_SA(indI))     
        OldTime=NewTime; OldK=indK;
        disp(' ')
    end
end
NewTime=toc;
       
% cut end
Ksat_MC(indK+1:end)=[]; c_MC(indK+1:end)=[]; tsub_MC(indK+1:end)=[]; z_MC(indK+1:end)=[]; 
EVAL(indK+1:end)=[]; T_SA_acc(indK+1:end)=[];

