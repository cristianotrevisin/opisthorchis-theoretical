% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;
% read OCN
set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')
OCN = build_OCN("OCN_A.mat",30*10000*10000);
% Add attributes


%figure; draw_OCN(OCN,NaN); hold on
%for sc = 1:OCN.nNodes; text(OCN.geometry.SCX(sc)/OCN.cellsize,OCN.geometry.SCY(sc)/OCN.cellsize,num2str(sc),'Color','k'); end;

par = common_parameters();

setup = build_setup(OCN,par,33*800*1000,'seed',3108);



outlet = find(OCN.SC_AccArea == max(OCN.SC_AccArea));
seeding_node = find(OCN.distW(outlet,:)==max(OCN.distW(outlet,:)));

y0 = zeros(OCN.nNodes,4);
y0(seeding_node,:) = [1/setup.H(seeding_node) 0 0 0];

Time = 1:20*365;

%Test without seasonality
setup.period = ones(length(Time),1);
yNS = model_ODE(Time,setup.par,setup,y0');
%Test with seasonality
setup.period = 0.5+0.5*sin((Time'-365.25)/365.25*2*pi);
yWS = model_ODE(Time,setup.par,setup,y0);

WHNS = yNS(:,1:4:end); WHWS = yWS(:,1:4:end);
EENS = yNS(:,2:4:end); EEWS = yWS(:,2:4:end);
YSNS = yNS(:,3:4:end); YSWS = yWS(:,3:4:end);
IFNS = yNS(:,4:4:end); IFWS = yWS(:,4:4:end);

setup.A(isnan(setup.A)) = 0; setup.S(isnan(setup.S))=0; setup.F(isnan(setup.F)) = 0;
WNS = WHNS*setup.H/sum(setup.H);
WWS = WHWS*setup.H/sum(setup.H);

ENS = EENS*setup.H/sum(setup.A);
EWS = EEWS*setup.H/sum(setup.A);

SNS = YSNS*setup.H/sum(setup.S);
SWS = YSWS*setup.H/sum(setup.S);

FNS = IFNS*setup.H/sum(setup.F);
FWS = IFWS*setup.H/sum(setup.F);

minWN = min(WHNS,[],2); minWW = min(WHWS,[],2); minWN(minWN < 0) = 1e-50;
maxWN = max(WHNS,[],2); maxWW = min(WHWS,[],2); minWW(minWN < 0) = 1e-50;
minEN = min(EENS,[],2); minEW = min(EEWS,[],2);
maxEN = max(EENS,[],2); maxEW = min(EEWS,[],2);
minSN = min(YSNS,[],2); minSW = min(YSWS,[],2);
maxSN = max(YSNS,[],2); maxSW = min(YSWS,[],2);
minFN = min(IFNS,[],2); minFW = min(IFWS,[],2);
maxFN = max(IFNS,[],2); maxFW = min(IFWS,[],2);
%%
close all
figure();
tiledlayout(2,2)
for k = 1:4
    nexttile
    if k == 1
        title('Worm burden in humans')
        hold on
        PLT2 = WNS;
        PLT4 = WWS;
        %filler(Time/365,maxWN',minWN',[0.16862745098039217 0.17647058823529413 0.25882352941176473],0.2)
        %filler(Time/365,maxWW',minWW',[0.8509803921568627 0.01568627450980392 0.1607843137254902],0.2)
    elseif k == 2
        title('Egg concentration')
        PLT1 = EENS;
        PLT2 = ENS;
        PLT3 = EEWS;
        PLT4 = EWS;
    elseif k == 3
        title('Prevalence of infected snails')
        PLT1 = YSNS;
        PLT2 = SNS;
        PLT3 = YSWS;
        PLT4 = SWS;
    elseif k == 4
        title('Cyst burden in fish')
        PLT1 = IFNS;
        PLT2 = FNS;
        PLT3 = IFWS;
        PLT4 = FWS;
    end
  
    hold on
    % for i = 1:OCN.nNodes
    %     a(i) = plot(Time/365,PLT1(:,i),'linewidth',0.25);
    %     a(i).Color = [0.16862745098039217 0.17647058823529413 0.25882352941176473 .1];
    % end
    plot(Time/365,PLT2,'color',[0.16862745098039217 0.17647058823529413 0.25882352941176473 1],'linewidth',1);
    % for i = 1:OCN.nNodes
    %     b(i) = plot(Time/365,PLT3(:,i),'linewidth',0.25);
    %     b(i).Color = [0.8509803921568627 0.01568627450980392 0.1607843137254902 .1];
    % end
    plot(Time/365,PLT4,'color',[0.8509803921568627 0.01568627450980392 0.1607843137254902 1],'linewidth',1);
    %set(gca,'Yscale','log')
    if k == 1
        legend('No seasonality','With seasonality','location','southeast')
    end
    xlabel('Years')
    if k == 1
        ylim([1e-20 Inf])
    elseif k == 2
        ylim([1e-25 Inf])
    elseif k == 3
        ylim([1e-30 Inf])
    elseif k == 4
        ylim([1e-25 Inf])
    end
end