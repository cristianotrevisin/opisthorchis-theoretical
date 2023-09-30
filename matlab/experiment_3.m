% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;
% read OCN

OCN = build_OCN("OCN_A.mat",30*10000*10000);
% Add attributes


% figure; draw_OCN(OCN,NaN); hold on
% for sc = 1:OCN.nNodes; text(OCN.geometry.SCX(sc)/OCN.cellsize,OCN.geometry.SCY(sc)/OCN.cellsize,num2str(sc),'Color','k'); end;

par = common_parameters();

setup = build_setup(OCN,par,33*800*1000,'seed',3108);
setup.T = eye(setup.nNodes);

outlet = find(OCN.SC_AccArea == max(OCN.SC_AccArea));
seeding_node = find(OCN.distW(outlet,:)==max(OCN.distW(outlet,:)));

y0 = zeros(OCN.nNodes,4);
y0(seeding_node,:) = [1/setup.H(seeding_node) 0 0 0];

Time = 1:100*365;

setup.period = ones(length(Time),1);

%% TEST DOWNSTREAM ONLY

setup.par.lambda_FD = 0.01;
par.lambda_FU = 0.00;

Times_Plot_Fig = [1 183 365 2*365 3*365 5*365 10*365 15*365 20*365 length(Time)];

Times_Plot_Map = 365*[182/365 1 2 4 7 10];
colorMap_IN = [linspace(016/256, 179/256,11); ...
    linspace(101/256, 021/256,11); ...
    linspace(171/256, 041/256,11)]';

colorMap_MP = [ones(256,1)';linspace(1,0,256);linspace(1,0,256)]';

y = model_ODE(Time,setup.par,setup,y0');
WH = y(:,1:4:end);
figure();plot(WH)


figure
for tt = 1:length(Times_Map)
        subplot(3,2,tt)
        draw_OCN(OCN,WH(Times_Plot_Map(tt),:)','Borders_Color','black')

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
