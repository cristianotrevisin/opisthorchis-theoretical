% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;
% read OCN

OCN = build_OCN("OCN_A.mat");
% Add attributes
OCN.thrA = 30*10000*10000; 

%figure; draw_OCN(OCN,NaN); hold on
%for sc = 1:OCN.nNodes; text(OCN.geometry.SCX(sc)/OCN.cellsize,OCN.geometry.SCY(sc)/OCN.cellsize,num2str(sc),'Color','k'); end;

par = common_parameters();
par.dF = 10;
par.dS = 30;
par.c = 3.0280e-08;
par.D = 10000;
par.lambda_FD = 0.01;
par.lambda_FU = 0.005;

[setup] = build_setup(OCN,par,33*800*1000,'seed',3108);


% Add parameter after volume
betaHS = 9.16e-11;
par.beta_E = betaHS*mean(setup.V)*par.mu_E/(par.rho_E-mean(setup.V)*betaHS*3102/7.5/0.5);


y0 = zeros(OCN.nNodes,5);
y0(:,4) = setup.KF;
y0(6,:) = [1 0 0 setup.KF(6) 0];

Time = 1:100*365;
setup.nNodes = OCN.nNodes;

%Test without seasonality
setup.period = ones(length(Time),1);
yNS = model_ODE(Time,par,setup,y0');
%Test with seasonality
setup.period = 0.625+0.375*sin((Time'-365.25)/365.25);
setup.period = 0.5+0.5*sin((Time'-365.25)/365.25);
yWS = model_ODE(Time,par,setup,y0');

WHNS = yNS(:,1:5:end);
WHWS = yWS(:,1:5:end);

WNS = WHNS*setup.H/sum(setup.H);
WWS = WHWS*setup.H/sum(setup.H);
%%
figure();
hold on
for i = 1:OCN.nNodes
    a(i) = plot(Time/365,WHNS(:,i),'linewidth',0.25);
    a(i).Color = [0.16862745098039217 0.17647058823529413 0.25882352941176473 .1];
end
plot(Time/365,WNS,'color',[0.16862745098039217 0.17647058823529413 0.25882352941176473 1],'linewidth',1);
for i = 1:OCN.nNodes
    b(i) = plot(Time/365,WHWS(:,i),'linewidth',0.25);
    b(i).Color = [0.8509803921568627 0.01568627450980392 0.1607843137254902 .1];
end
plot(Time/365,WWS,'color',[0.8509803921568627 0.01568627450980392 0.1607843137254902 1],'linewidth',1);
set(gca,'Yscale','log')
legend('','No seasonality','','With seasonality')
ylim([1e-8 1e8])
ylabel('I^H')
xlabel('Years')