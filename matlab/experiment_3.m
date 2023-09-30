% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;
% read OCN

OCN = build_OCN("OCN_A.mat");
% Add attributes
OCN.thrA = 30*10000*10000; 

figure; draw_OCN(OCN,NaN); hold on
for sc = 1:OCN.nNodes; text(OCN.geometry.SCX(sc)/OCN.cellsize,OCN.geometry.SCY(sc)/OCN.cellsize,num2str(sc),'Color','k'); end;

par = common_parameters();


[setup] = build_setup(OCN,par,33*800*1000,'seed',3108);
setup.T = eye(setup.nNodes);

% Add parameter after volume
betaHS = 9.16e-11;
par.beta_E = betaHS*mean(setup.V)*par.mu_E/(par.rho_E-mean(setup.V)*betaHS*3102/7.5/0.5);

F = find_fish_equilibrium(setup.KF,setup.W,par,setup.chi);

y0 = zeros(OCN.nNodes,5);
y0(:,4) = F;
y0(6,:) = [1 0 0 setup.KF(6) 0];

Time = 1:100*365;
setup.nNodes = OCN.nNodes;


%% TEST DOWNSTREAM ONLY

par.lambda_FD = 0.01;
par.lambda_FU = 0.005;

Times = [365:365:365*10];

Times_Map = 365*[182/365 1 2 4 7 10];
colorMap_IN = [linspace(016/256, 179/256,11); ...
    linspace(101/256, 021/256,11); ...
    linspace(171/256, 041/256,11)]';

colorMap_MP = [ones(256,1)';linspace(1,0,256);linspace(1,0,256)]';

y = model_ODE(Time,par,setup,y0');
WH = y(:,1:5:end);

%WH(:,6) = NaN;

figure
for tt = 1:length(Times_Map)
        subplot(3,2,tt)
        draw_OCN(OCN,WH(Times(tt+1),:)','Borders_Color','black')

        set(gca,'ColorScale','log')
        colorbar
        clim([1e-5 5e6]);
        colormap(colorMap_MP)
        colorbar( 'off' ) 
end

dist = OCN.distW(6,:);

MTR = dist;
for i = 1:length(Times)
    ind = Time==Times(i);
    MTR = [MTR;WH(ind,:)];
end
MTR = MTR';


idx_delete = find(MTR(:,end)==0);
MTR(idx_delete,:) = [];

MTR = sortrows(MTR);



figure();
hold on
for tt = 1:length(Times)
    plot(MTR(:,1)/1000,MTR(:,tt+1),'color',colorMap_IN(tt,:),...
        'linestyle','-','marker','*','linewidth',1)
end
set(gca,'Yscale','log')
xlabel('Distance to first infected node [km]')
ylabel('I^H')
legend(num2str(Times'/365))
