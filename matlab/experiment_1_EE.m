% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;
% read OCN

OCN = build_OCN("OCN_A.mat");
% Add attributes
OCN.thrA = 30*10000*10000; 


par = common_parameters();
par.dF = 10;
par.dS = 30;
par.D = 100000;
par.lambda_FU=0;
par.lambda_FD=0;
par.lambda_ED=0;


[setup_RA,par_RA] = build_setup(OCN,par,33*800*1000,'seed',3108);
[setup_DA,par_DA] = build_setup(OCN,par,33*800*1000,'seed',3108,'DownstreamAccumulation',true);
[setup_UNI,parU] = build_setup(OCN,par,33*800*1000,'seed',3108,'Unify',true);

% Prevent transportation of raw fish outside community
setup_RA.T = eye(OCN.nNodes); %setup.T(eye(OCN.nNodes)==1) = 1-setup.sigma;
setup_DA.T = eye(OCN.nNodes); %setup.T(eye(OCN.nNodes)==1) = 1-setup.sigma;



% Now check endemic equilibrium in nodes (random allocation)
eqi_RA = zeros(OCN.nNodes,5);
for nn = 1:OCN.nNodes
    tmp = max(find_EE(par_RA,par_RA.c*setup_RA.H(nn),setup_RA.H(nn),...
        setup_RA.S(nn),setup_RA.KF(nn),setup_RA.V(nn)));
    if  isempty(tmp) || tmp(1)==0
        eqi_RA(nn,:)=nan;
    else
        eqi_RA(nn,:)=tmp;
    end
end
EE_reached_RA = ~isnan(eqi_RA(:,1));
sum(EE_reached_RA)


% Now check endemic equilibrium in nodes (downstream accumulation)
eqi_DA = zeros(OCN.nNodes,5);
for nn = 1:OCN.nNodes
    tmp = max(find_EE(par_DA,par_DA.c*setup_DA.H(nn),setup_DA.H(nn),...
        setup_DA.S(nn),setup_DA.KF(nn),setup_DA.V(nn)));
    if  isempty(tmp) || tmp(1)==0
        eqi_DA(nn,:)=nan;
    else
        eqi_DA(nn,:)=tmp;
    end
end
EE_reached_DA = ~isnan(eqi_DA(:,1));
sum(EE_reached_DA)

% Check endemic equilibrium in whole system
eqi = find_EE(parU,parU.c*setup_UNI.H,setup_UNI.H,setup_UNI.S,setup_UNI.KF,setup_UNI.V);

% Running simulations

% Random allocation
y0_RA = zeros(5,OCN.nNodes);
y0_RA(4,:) = setup_RA.KF';
y0_RA(1,:) = 1./setup_RA.H';

% Downstream accumulation
y0_DA = zeros(5,OCN.nNodes);
y0_DA(4,:) = setup_DA.KF';
y0_DA(1,:) = 1./setup_DA.H';

% Unified
y0_UNI = zeros(5,1);
y0_UNI(4) = setup_UNI.KF;
y0_UNI(1) = 1./setup_UNI.H;


% Launch simulations
Time = 1:1000*365;
y_RA = model_ODE(Time,par_RA,setup_RA,y0_RA);
y_DA = model_ODE(Time,par_DA,setup_DA,y0_DA);
y_UNI = model_ODE(Time,parU,setup_UNI,y0_UNI);


% Agglomerate results
WH_SE_RA = y_RA(:,1:5:end)*setup_RA.H/(sum(setup_RA.H));
E_SE_RA = y_RA(:,2:5:end)*setup_RA.V/(sum(setup_RA.V));
YS_SE_RA = y_RA(:,3:5:end)*setup_RA.S/(sum(setup_RA.S));
F_SE_RA = sum(y_RA(:,4:5:end),2);
IH_SE_RA = sum(y_RA(:,4:5:end).*y_RA(:,5:5:end),2)./F_SE_RA;

WH_RA = y_RA(:,1:5:end);
WH_DA = y_DA(:,1:5:end);

WH_SE_DA = y_DA(:,1:5:end)*setup_DA.H/(sum(setup_DA.H));
E_SE_DA = y_DA(:,2:5:end)*setup_DA.V/(sum(setup_DA.V));
YS_SE_DA = y_DA(:,3:5:end)*setup_DA.S/(sum(setup_DA.S));
F_SE_DA = sum(y_DA(:,4:5:end),2);
IH_SE_DA = sum(y_DA(:,4:5:end).*y_DA(:,5:5:end),2)./F_SE_DA;



figure
loglog(Time/365,WH_RA(:,EE_reached_RA==1),'color','#cc2936');
hold on
loglog(Time/365,WH_RA(:,EE_reached_RA==0),'color','#08415c');
xlim([1/10 Time(end)/365])
ylim([1/max(setup_RA.H) max(WH_RA,[],'all')])
ylabel('I^H')
xlabel('Years')


figure
draw_OCN(OCN,EE_reached_RA,'binary',true)
% for sc = 1:OCN.nNodes
%     text(OCN.geometry.SCX(sc)/OCN.cellsize,OCN.geometry.SCY(sc)/OCN.cellsize,num2str(sc),'Color','white')
% end


figure
loglog(Time/365,WH_DA(:,EE_reached_DA==1),'color','#cc2936');
hold on
loglog(Time/365,WH_DA(:,EE_reached_DA==0),'color','#08415c');
xlim([1/10 Time(end)/365])
ylim([1/max(setup_DA.H) max(WH_DA,[],'all')])
ylabel('I^H')
xlabel('Years')

figure
draw_OCN(OCN,EE_reached_DA,'binary',true)
% for sc = 1:OCN.nNodes
%     text(OCN.geometry.SCX(sc)/OCN.cellsize,OCN.geometry.SCY(sc)/OCN.cellsize,num2str(sc),'Color','white')
% end


%%
figure()

tiledlayout(2,2)

nexttile()
plot(Time/365,y_UNI(:,1))
hold on
plot(Time/365,WH_SE_RA)
plot(Time/365,WH_SE_DA)
plot(Time/365,ones(length(Time),1)*max(eqi(:,1)),'--r')
legend('At-large','RA','DA','EE At-large','Orientation','horizontal')
ylabel('I^H [w/h]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
set(gca,'XScale','log','YScale','log')
box off
grid minor 

% subplot(2,2,2)
% plot(Time/365,y_UNI(:,2))
% hold on
% plot(Time/365,E_SE_RA)
% plot(Time/365,E_SE_DA)
% plot(Time/365,ones(length(Time),1)*max(eqi(:,2)),'--r')
% ylabel('E [e/m^3]')
% xlim([Time(1)/365 Time(end)/365])
% set(gca,'XTickLabel',[])
% set(gca,'XScale','log','YScale','log')
% box off

nexttile()
plot(Time/365,100*y_UNI(:,3))
hold on
plot(Time/365,YS_SE_RA)
plot(Time/365,YS_SE_DA)
plot(Time/365,100*ones(length(Time),1)*max(eqi(:,3)),'--r')
ylabel('Y^S [%]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
set(gca,'XScale','log','YScale','log')
box off
grid minor 


nexttile()
plot(Time/365,y_UNI(:,4))
hold on
plot(Time/365,F_SE_RA)
plot(Time/365,F_SE_DA)
plot(Time/365,ones(length(Time),1)*max(eqi(:,4)),'--r')
ylabel('F [-]')
xlim([Time(1)/365 Time(end)/365])
%set(gca,'XTickLabel',[])
set(gca,'XScale','log','YScale','log')
box off
xlabel('Years')
grid minor

nexttile()
plot(Time/365,y_UNI(:,5))
hold on
plot(Time/365,IH_SE_RA)
plot(Time/365,IH_SE_DA)
plot(Time/365,ones(length(Time),1)*max(eqi(:,5)),'--r')
ylabel('I^F [m/f]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XScale','log','YScale','log')
xlabel('Years')
box off
grid minor