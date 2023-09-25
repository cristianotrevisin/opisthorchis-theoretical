% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;
% read OCN

OCN_A = build_OCN("OCN_A.mat");
OCN_B = build_OCN("OCN_B.mat");
OCN_C = build_OCN("OCN_C.mat");
% Add attributes
OCN_A.thrA = 30*10000*10000; 
OCN_B.thrA = 30*10000*10000; 
OCN_C.thrA = 30*10000*10000; 

par = common_parameters();
par.dF = 10;
par.dS = 30;
par.lambda_FU=0;
par.lambda_FD=0;
par.lambda_ED=0;
par.D = 100000;


[setup,par] = build_setup(OCN_A,par,33*800*1000,'seed',3108);

[setup,par] = build_setup(OCN_A,par,33*800*1000,'seed',3108,'Unify',true);

% figure(); imagesc(setup.T-diag(diag(setup.T)))


figure();

Time = 1:20*365;


y0 = [1/setup.H; 0; 0; setup.KF; 0];


eqi = find_EE(par,par.c*setup.H,setup.H,setup.S,setup.KF,setup.V);

y = model_ODE(Time,par,setup,y0);

figure()
subplot(5,1,1)
plot(Time/365,y(:,1))
hold on
plot(Time/365,ones(length(Time),1)*max(eqi(:,1)),'--r')
ylabel('I^H [w/h]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
box off

subplot(5,1,2)
plot(Time/365,y(:,2))
hold on
plot(Time/365,ones(length(Time),1)*max(eqi(:,2)),'--r')
ylabel('E [e/m^3]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
box off

subplot(5,1,3)
plot(Time/365,100*y(:,3))
hold on
plot(Time/365,100*ones(length(Time),1)*max(eqi(:,3)),'--r')
ylabel('Y^S [%]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
box off

subplot(5,1,4)
plot(Time/365,y(:,4))
hold on
plot(Time/365,ones(length(Time),1)*max(eqi(:,4)),'--r')
ylabel('F [-]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
box off

subplot(5,1,5)
plot(Time/365,y(:,5))
hold on
plot(Time/365,ones(length(Time),1)*max(eqi(:,5)),'--r')
ylabel('I^F [m/f]')
xlim([Time(1)/365 Time(end)/365])
xlabel('Years')
box off

%% Now check endemic equilibrium in nodes
eqi_NN = zeros(OCN_A.nNodes,5);

for nn = 1:OCN_A.nNodes
    tmp = max(find_EE(par,par.c*setup_N.H(nn),setup_N.H(nn),...
        setup_N.S(nn),setup_N.KF(nn),setup_N.V(nn)));
    if isempty(tmp) || tmp(1)==0
        eqi_NN(nn,:)=nan;
    else
        eqi_NN(nn,:)=tmp;
    end
end

figure()
draw_OCN(OCN_A,eqi_NN(:,1));
colorbar
for sc = 1:OCN_A.nNodes
    plot(OCN_A.geometry.SCX(sc)/OCN_A.cellsize,OCN_A.geometry.SCY(sc)/OCN_A.cellsize,'.r','MarkerSize',0.5+1.5*log(setup.H(sc)))
    text(OCN_A.geometry.SCX(sc)/OCN_A.cellsize,OCN_A.geometry.SCY(sc)/OCN_A.cellsize,num2str(sc),'Color','k')
end

setup_N.T = eye(OCN_A.nNodes); setup_N.T(eye(OCN_A.nNodes)==1) = 1-setup_N.sigma;
y0 = zeros(OCN_A.nNodes,5);
y0(:,4) = setup.KF;
y0(:,1) = 1./setup_N.H;

Time = 1:220*365;


par.lambda_FU=0;
par.lambda_FD=0;
par.lambda_ED=0;
setup_N.nNodes = OCN_A.nNodes;
y = model_ODE(Time,par,setup_N,y0');

WH = y(:,1:5:end);

figure();
semilogy(Time/365,WH')
for nn = 1:OCN_A.nNodes
    text(Time(end)/365,WH(end,nn),num2str(nn));
end
%%
% Now check what happens with downstream ACC
eqi_ND = zeros(OCN_A.nNodes,5);

[setup_N,par] = build_setup(OCN_A,par,33*800*1000,'DownstreamAccumulation', true,'seed',3108);
for nn = 1:OCN_A.nNodes
    tmp = max(find_EE(par,par.c*setup_N.H(nn),setup_N.H(nn),...
        setup_N.S(nn),setup_N.KF(nn),setup_N.V(nn)));
    if tmp(1)==0 || isempty(tmp)
        eqi_ND(nn,:)=nan;
    else
        eqi_ND(nn,:)=tmp;
    end
end

figure()
draw_OCN(OCN_A,eqi_ND(:,1));
colorbar
%% SENSITIVITY ANALYSIS

[ups, downs, labels] = sensitivity_analysis_EE(OCN_A,par,...
    {'\lambda_FU';'\lambda_FD';'\lambda_ED';'D'});

figure
barh(1:length(labels),downs(1,:)*100)
hold on
barh(1:length(labels),ups(1,:)*100)
xlabel('[%] prevalence variation')
box off
ax = gca;
%set(h,'ycolor','w')
ax.YColor = 'w';
ax.YAxis.Label.Color='k';
yticks(1:length(labels))
yticklabels(labels)
set(gca,'YDir','reverse')
