% In here we want to investigate the existance of a local equilibrium


clc; close all; clearvars;

par = common_parameters();

Time = 1:100*365; 

OCN = build_OCN("OCN_A.mat",30*10000*10000);
sd = 3108;


% Load setups

setup = build_setup(OCN,par,33*800*1000,'seed',sd);


% Prevent transportation of raw fish outside community
setup.T(eye(OCN.nNodes)==1) = 1-setup.sigma;


% Now check endemic equilibrium in nodes (random allocation)
eqi_RA = zeros(OCN.nNodes,4);
for nn = 1:OCN.nNodes
    tmp = max(find_EE(setup.par,setup.par.c*setup.H(nn),setup.H(nn),...
        setup.S(nn),setup.F(nn),setup.A(nn),setup.par.theta_C));
    if  isempty(tmp) || tmp(1)==0
        eqi_RA(nn,:)=0;
    else
        eqi_RA(nn,:)=tmp;
    end
end
EE_reached_RA = eqi_RA(:,1)>0;
sum(EE_reached_RA)



% Random allocation
y0_RA = zeros(4,OCN.nNodes);
y0_RA(1,:) = 1./setup.H';

% Launch simulations
setup.period = ones(length(Time),1); 
y_RA = model_ODE(Time,setup.par,setup,y0_RA);