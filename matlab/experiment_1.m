%PRELIMINARY COMPARISON

clc; close all; clearvars;

set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

par = common_parameters();

Time = 1:100*365; 

plt_RA = zeros(length(Time),3);
plt_DA = zeros(length(Time),3);
plt_UNI = zeros(length(Time),3);

% USER SHOULD SELECT WHICH OCNMAP TO USE
ocnmap = 3;

if ocnmap == 1
    OCN = build_OCN("OCN_A.mat",30*10000*10000);
    sd = 3108;
elseif ocnmap == 2
    OCN = build_OCN("OCN_B.mat",30*10000*10000);
    sd = 2507;
elseif ocnmap == 3
    OCN = build_OCN("OCN_C.mat",30*10000*10000);
    sd = 0808;
end

% Load setups
par.lambda_F = 0;
setup_RA = build_setup(OCN,par,33*800*1000,'seed',sd);
setup_DA = build_setup(OCN,par,33*800*1000,'seed',sd,'DownstreamAccumulation',true);


temp_T_RA = setup_RA.T; temp_T_DA = setup_DA.T;

% We want a vector, not a matrix, for compute_EE
setup_RA.T = (min(setup_RA.par.U*setup_RA.H,setup_RA.chi.*setup_RA.F));
setup_DA.T = (min(setup_DA.par.U*setup_DA.H,setup_DA.chi.*setup_DA.F));


% Now check endemic equilibrium in nodes (random allocation)
eqi_RA = compute_EE(setup_RA);
EE_reached_RA = eqi_RA(:,1)>0;
sum(EE_reached_RA)
    

% Now check endemic equilibrium in nodes (downstream accumulation)
eqi_DA = compute_EE(setup_DA);
EE_reached_DA = eqi_DA(:,1)>0;
sum(EE_reached_DA)

  
% Running simulations
setup_RA.A(isnan(setup_RA.A)) = 0; setup_DA.A(isnan(setup_DA.A)) = 0;
% Random allocation
y0_RA = zeros(4,OCN.nNodes);
y0_RA(1,:) = 0.1;

% Downstream accumulation
y0_DA = zeros(4,OCN.nNodes);
y0_DA(1,:) = 0.1;
    
% Restore diagonal matrix
setup_RA.T = diag(min(setup_RA.par.U*setup_RA.H,setup_RA.chi.*setup_RA.F));
setup_DA.T = diag(min(setup_DA.par.U*setup_DA.H,setup_DA.chi.*setup_DA.F));

% Launch simulations

% Without mobility
y_RA_1 = model_ODE(Time,setup_RA.par,setup_RA,y0_RA);
y_DA_1 = model_ODE(Time,setup_DA.par,setup_DA,y0_DA);

% Only fish connectivity
setup_RA.par.lambda_F = 0.01; setup_DA.par.lambda_F = 0.01;
setup_RA.W = calculateWRecursive(1, setup_RA.X+setup_RA.X', setup_RA.F, zeros(size(setup_RA.X)), false(setup_RA.nNodes,1),1);   
setup_RA.W = setup_RA.W*setup_RA.par.lambda_F/mean(setup_RA.W(:)>0,'all');

setup_DA.W = calculateWRecursive(1, setup_DA.X+setup_DA.X', setup_DA.F, zeros(size(setup_DA.X)), false(setup_DA.nNodes,1),1);   
setup_DA.W = setup_RA.W*setup_DA.par.lambda_F/mean(setup_DA.W(:)>0,'all');


y_RA_2 = model_ODE(Time,setup_RA.par,setup_RA,y0_RA);
y_DA_2 = model_ODE(Time,setup_DA.par,setup_DA,y0_DA);

% Only fish market
setup_RA.par.lambda_F = 0; setup_DA.par.lambda_F = 0;
setup_RA.T = temp_T_RA; setup_DA.T = temp_T_DA; 
setup_RA.W = calculateWRecursive(1, setup_RA.X+setup_RA.X', setup_RA.F, zeros(size(setup_RA.X)), false(setup_RA.nNodes,1),1);   
setup_RA.W = setup_RA.W*setup_RA.par.lambda_F/mean(setup_RA.W(:)>0,'all');

setup_DA.W = calculateWRecursive(1, setup_DA.X+setup_DA.X', setup_DA.F, zeros(size(setup_DA.X)), false(setup_DA.nNodes,1),1);   
setup_DA.W = setup_RA.W*setup_DA.par.lambda_F/mean(setup_DA.W(:)>0,'all');

y_RA_3 = model_ODE(Time,setup_RA.par,setup_RA,y0_RA);
y_DA_3 = model_ODE(Time,setup_DA.par,setup_DA,y0_DA);

% Both kinds of mobility
setup_RA.par.lambda_F = 0.01; setup_DA.par.lambda_F = 0.01;

setup_RA.W = calculateWRecursive(1, setup_RA.X+setup_RA.X', setup_RA.F, zeros(size(setup_RA.X)), false(setup_RA.nNodes,1),1);   
setup_RA.W = setup_RA.W*setup_RA.par.lambda_F/mean(setup_RA.W(:)>0,'all');

setup_DA.W = calculateWRecursive(1, setup_DA.X+setup_DA.X', setup_DA.F, zeros(size(setup_DA.X)), false(setup_DA.nNodes,1),1);   
setup_DA.W = setup_RA.W*setup_DA.par.lambda_F/mean(setup_DA.W(:)>0,'all');

setup_RA.T = temp_T_RA; setup_DA.T = temp_T_DA; 
y_RA_4 = model_ODE(Time,setup_RA.par,setup_RA,y0_RA);
y_DA_4 = model_ODE(Time,setup_DA.par,setup_DA,y0_DA);
    


% Agglomerate results
plt_RA(:,1,1) = y_RA_1(:,1:3:end)*setup_RA.H/(sum(setup_RA.H));
plt_RA(:,2,1) = y_RA_1(:,2:3:end)*setup_RA.S/(sum(setup_RA.S));
plt_RA(:,3,1) = y_RA_1(:,3:3:end)*setup_RA.F/(sum(setup_RA.F));

plt_RA(:,1,2) = y_RA_2(:,1:3:end)*setup_RA.H/(sum(setup_RA.H));
plt_RA(:,2,2) = y_RA_2(:,2:3:end)*setup_RA.S/(sum(setup_RA.S));
plt_RA(:,3,2) = y_RA_2(:,3:3:end)*setup_RA.F/(sum(setup_RA.F));

plt_RA(:,1,3) = y_RA_3(:,1:3:end)*setup_RA.H/(sum(setup_RA.H));
plt_RA(:,2,3) = y_RA_3(:,2:3:end)*setup_RA.S/(sum(setup_RA.S));
plt_RA(:,3,3) = y_RA_3(:,3:3:end)*setup_RA.F/(sum(setup_RA.F));

plt_RA(:,1,4) = y_RA_4(:,1:3:end)*setup_RA.H/(sum(setup_RA.H));
plt_RA(:,2,4) = y_RA_4(:,2:3:end)*setup_RA.S/(sum(setup_RA.S));
plt_RA(:,3,4) = y_RA_4(:,3:3:end)*setup_RA.F/(sum(setup_RA.F));


plt_DA(:,1,1) = y_DA_1(:,1:3:end)*setup_DA.H/(sum(setup_DA.H));
plt_DA(:,2,1) = y_DA_1(:,2:3:end)*setup_DA.S/(sum(setup_DA.S));
plt_DA(:,3,1) = y_DA_1(:,3:3:end)*setup_DA.F/(sum(setup_DA.F));

plt_DA(:,1,2) = y_DA_2(:,1:3:end)*setup_DA.H/(sum(setup_DA.H));
plt_DA(:,2,2) = y_DA_2(:,2:3:end)*setup_DA.S/(sum(setup_DA.S));
plt_DA(:,3,2) = y_DA_2(:,3:3:end)*setup_DA.F/(sum(setup_DA.F));

plt_DA(:,1,3) = y_DA_3(:,1:3:end)*setup_DA.H/(sum(setup_DA.H));
plt_DA(:,2,3) = y_DA_3(:,2:3:end)*setup_DA.S/(sum(setup_DA.S));
plt_DA(:,3,3) = y_DA_3(:,3:3:end)*setup_DA.F/(sum(setup_DA.F));

plt_DA(:,1,4) = y_DA_4(:,1:3:end)*setup_DA.H/(sum(setup_DA.H));
plt_DA(:,2,4) = y_DA_4(:,2:3:end)*setup_DA.S/(sum(setup_DA.S));
plt_DA(:,3,4) = y_DA_4(:,3:3:end)*setup_DA.F/(sum(setup_DA.F));

WH_RA_1 = y_RA_1(:,1:3:end);
WH_DA_1 = y_DA_1(:,1:3:end);
WH_RA_2 = y_RA_2(:,1:3:end);
WH_DA_2 = y_DA_2(:,1:3:end);
WH_RA_3 = y_RA_3(:,1:3:end);
WH_DA_3 = y_DA_3(:,1:3:end);
WH_RA_4 = y_RA_4(:,1:3:end);
WH_DA_4 = y_DA_4(:,1:3:end);


%% Generate figures
figure()

colors = ["#353535", "#edae49", "#d1495b", "#00798c"];
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile([1 2])
hold on
for scenario = 1:4
    plot(Time/365,plt_RA(:,1,scenario),'color',colors(scenario),'LineStyle','-','LineWidth',1)
    plot(Time/365,plt_DA(:,1,scenario),'color',colors(scenario),'LineStyle','--','LineWidth',1)
end
ylabel('I^H [worms/person]')
%set(gca,'XTick',[0.01 0.1 1 10 100])
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
%set(gca,'XScale','log')
%legend('0', "", 'Fish mobility',"",'Fish market',"",'Both')
box off
set(gca,'FontSize',9)


nexttile([1 2])
hold on
for scenario = 1:4
    plot(Time/365,plt_RA(:,2,scenario),'color',colors(scenario),'LineStyle','-','LineWidth',1)
    plot(Time/365,plt_DA(:,2,scenario),'color',colors(scenario),'LineStyle','--','LineWidth',1)
end
ylabel('Y^S')
%set(gca,'XTick',[0.01 0.1 1 10 100])
xlim([Time(1)/365 Time(end)/365])
%set(gca,'XScale','log')
box off
xlabel('Years')
set(gca,'FontSize',9)

nexttile([1 2])
hold on
for scenario = 1:4
    plot(Time/365,plt_RA(:,3,scenario),'color',colors(scenario),'LineStyle','-','LineWidth',1)
    plot(Time/365,plt_DA(:,3,scenario),'color',colors(scenario),'LineStyle','--','LineWidth',1)
end
ylabel('I^F [cysts/fish]')
%set(gca,'XTick',[0.01 0.1 1 10 100])
xlim([Time(1)/365 Time(end)/365])
%set(gca,'XScale','log')
xlabel('Years')
box off
set(gca,'FontSize',9)

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [17 12]);
set(findall(gcf,'-property','FontSize'),'FontSize',9)

figure
draw_OCN(OCN,EE_reached_RA,'binary',true)
% 
figure
draw_OCN(OCN,EE_reached_DA,'binary',true)