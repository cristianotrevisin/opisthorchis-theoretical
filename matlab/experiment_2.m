% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;

set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

par = common_parameters();
Time = 1:25*365;
Ds = [1e4 5e4 1e5 5e5];

for ocnmap = 1:3
    if ocnmap == 1
        OCN = build_OCN("OCN_A.mat",30*10000*10000);
        sd = 3108;
        WH1 = zeros(length(Time),OCN.nNodes,length(Ds));
    elseif ocnmap == 2
        OCN = build_OCN("OCN_B.mat",30*10000*10000);
        sd = 2507;
        WH2 = zeros(length(Time),OCN.nNodes,length(Ds));
    elseif ocnmap == 3
        OCN = build_OCN("OCN_C.mat",30*10000*10000);
        sd = 0808;
        WH3 = zeros(length(Time),OCN.nNodes,length(Ds));
    end

   
    for scd = 1:length(Ds)
        par.D = Ds(scd);
        setup = build_setup(OCN,par,33*800*1000,'seed',sd);
        setup.period = ones(length(Time),1);
        surplus = setup.sigma.*setup.F;
    
        y0 = zeros(OCN.nNodes,4);
        y0(surplus==max(surplus),1) = 1/setup.H(surplus==max(surplus));
    
        y = model_ODE(Time,setup.par,setup,y0');

        if ocnmap == 1
            DIST1 = OCN.Dist(surplus==max(surplus),:);
    
            WH1(:,:,scd) = y(:,1:4:end).*setup.H';
        
        elseif ocnmap == 2
            DIST2 = OCN.Dist(surplus==max(surplus),:);

            WH2(:,:,scd) = y(:,1:4:end).*setup.H';

        elseif ocnmap == 3
            DIST3 = OCN.Dist(surplus==max(surplus),:);
    
            WH3(:,:,scd) = y(:,1:4:end).*setup.H';
        end
    end
end

%%
% PLOT
% Mix distances
DIST = [DIST1 DIST2 DIST3];

% Choose times
Times_Plot_Map = [30 180 5*365];
Times_Plot_Fig = [1 183 365 2*365 3*365 5*365 10*365 15*365 20*365 length(Time)];
colorMap_IN = [linspace(016/256, 179/256,11); ...
    linspace(101/256, 021/256,11); ...
    linspace(171/256, 041/256,11)]';



figure

for scd = 1:length(Ds)
    WH = [DIST;
        [squeeze(WH1(Times_Plot_Fig,:,scd))';
         squeeze(WH2(Times_Plot_Fig,:,scd))';
         squeeze(WH3(Times_Plot_Fig,:,scd))']'];

    WH(:,abs(WH(end,:))< 1e-10) = [];

    WH = sortrows(WH');
    nexttile()

    hold on
    for tt = 1:length(Times_Plot_Fig)
        plot(WH(:,1)/1000,exp(smoothdata(log10(WH(:,tt+1)),'smoothingfactor',1)),'color',colorMap_IN(tt,:),...
            'linestyle','-','linewidth',1)
        %plot(WH(:,1)/1000,exp(log10(WH(:,tt+1))),'color',colorMap_IN(tt,:),...
        %    'linestyle','-','linewidth',1)
    end
    set(gca,'Yscale','log')
    xlabel('Distance to first infected node [km]')
    ylabel('I^H')
    if scd == length(Ds)
        legend(num2str(Times_Plot_Fig'/365),'NumColumns',2)
    end
    ylim([1e-8 1e5])
    xlim([0 1000])
end


%%
OCN = build_OCN("OCN_A.mat",30*10000*10000);
        sd = 3108;
colorMap_MP = [ones(256,1)';linspace(1,0,256);linspace(1,0,256)]';
figure
tiledlayout(3,3)
cntr = 0;
for d = 1:3
    for tt = 1:length(Times_Plot_Map)
        nexttile()
        draw_OCN(OCN,WH1(Times_Plot_Map(tt),:,d)','Borders_Color','black')
        set(gca,'ColorScale','log')
        colorbar
        clim([1e-5 5e6]);
        colormap(colorMap_MP)
        colorbar( 'off' ) 
    end
end
