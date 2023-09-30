% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;

par = common_parameters();
Time = 1:25*365;
Ds = [1e4 5e4 1e5 5e5];

for ocnmap = 1:3
    if ocnmap == 1
        OCN = build_OCN("OCN_A.mat",30*10000*10000);
        sd = 3108;
        WH1_RA = zeros(length(Time),OCN.nNodes,length(Ds));
        WH1_DA = zeros(length(Time),OCN.nNodes,length(Ds));
    elseif ocnmap == 2
        OCN = build_OCN("OCN_B.mat",30*10000*10000);
        sd = 2507;
        WH2_RA = zeros(length(Time),OCN.nNodes,length(Ds));
        WH2_DA = zeros(length(Time),OCN.nNodes,length(Ds));
    elseif ocnmap == 3
        OCN = build_OCN("OCN_C.mat",30*10000*10000);
        sd = 0808;
        WH3_RA = zeros(length(Time),OCN.nNodes,length(Ds));
        WH3_DA = zeros(length(Time),OCN.nNodes,length(Ds));
    end

   
    for scd = 1:length(Ds)
        par.D = Ds(scd);
        setup_RA = build_setup(OCN,par,33*800*1000,'seed',sd);
        setup_DA = build_setup(OCN,par,33*800*1000,'seed',sd,'DownstreamAccumulation',true);
        setup_RA.period = ones(length(Time),1); setup_DA.period = ones(length(Time),1);
        surplus_RA = setup_RA.sigma.*setup_RA.F;
        surplus_DA = setup_DA.sigma.*setup_DA.F;
    
        y0_RA = zeros(OCN.nNodes,4);
        y0_RA(surplus_RA==max(surplus_RA),1) = 1/setup_RA.H(surplus_RA==max(surplus_RA));
    
        y0_DA = zeros(OCN.nNodes,4);
        y0_DA(surplus_RA==max(surplus_DA),1) = 1/setup_DA.H(surplus_DA==max(surplus_DA));

        y_RA = model_ODE(Time,setup_RA.par,setup_RA,y0_RA');
        y_DA = model_ODE(Time,setup_DA.par,setup_DA,y0_DA');

        if ocnmap == 1
            DIST1_RA = OCN.Dist(surplus_RA==max(surplus_RA),:);
            DIST1_DA = OCN.Dist(surplus_DA==max(surplus_DA),:);
    
            WH1_RA(:,:,scd) = y_RA(:,1:4:end);
            WH1_DA(:,:,scd) = y_DA(:,1:4:end);
        
        elseif ocnmap == 2
            DIST2_RA = OCN.Dist(surplus_RA==max(surplus_RA),:);
            DIST2_DA = OCN.Dist(surplus_DA==max(surplus_DA),:);
    
            WH2_RA(:,:,scd) = y_RA(:,1:4:end);
            WH2_DA(:,:,scd) = y_DA(:,1:4:end);

        elseif ocnmap == 3
            DIST3_RA = OCN.Dist(surplus_RA==max(surplus_RA),:);
            DIST3_DA = OCN.Dist(surplus_DA==max(surplus_DA),:);
    
            WH3_RA(:,:,scd) = y_RA(:,1:4:end);
            WH3_DA(:,:,scd) = y_DA(:,1:4:end);
        end
    end
end


% PLOT
% Mix distances
DIST_RA = [DIST1_RA DIST2_RA DIST3_RA];
DIST_DA = [DIST1_DA DIST2_DA DIST3_DA];




% Choose times
Times_Plot_Map = [30 365 20*365];
Times_Plot_Fig = [1 183 365 2*365 3*365 5*365 10*365 15*365 20*365 length(Time)];
colorMap_IN = [linspace(016/256, 179/256,11); ...
    linspace(101/256, 021/256,11); ...
    linspace(171/256, 041/256,11)]';




%%
figure

for scd = 1:length(Ds)
    WH = [DIST_RA DIST_DA;
        [squeeze(WH1_RA(Times_Plot_Fig,:,scd))';
         squeeze(WH2_RA(Times_Plot_Fig,:,scd))';
         squeeze(WH3_RA(Times_Plot_Fig,:,scd))';
         squeeze(WH1_DA(Times_Plot_Fig,:,scd))';
         squeeze(WH2_DA(Times_Plot_Fig,:,scd))';
         squeeze(WH3_DA(Times_Plot_Fig,:,scd))']'];

    WH(:,abs(WH(end,:))< 1e-10) = [];

    WH = sortrows(WH');
    nexttile()

    hold on
    for tt = 1:length(Times_Plot_Fig)
        plot(WH(:,1)/1000,exp(smoothdata(log10(WH(:,tt+1)),'smoothingfactor',1)),'color',colorMap_IN(tt,:),...
            'linestyle','-','linewidth',1)
    end
    set(gca,'Yscale','log')
    xlabel('Distance to first infected node [km]')
    ylabel('I^H')
    %legend(num2str(Times_Plot_Fig'/365))
    
end


%%

colorMap_MP = [ones(256,1)';linspace(1,0,256);linspace(1,0,256)]';
figure
tiledlayout(2,3)
cntr = 0;
for d = Ds(2):Ds(3)
    for tt = 1:length(Times_Plot_Map)
        cntr = cntr+1;
        subplot(4,4,cntr)
        draw_OCN(OCN,WH1_DA(tt,:,d)','Borders_Color','black')
        set(gca,'ColorScale','log')
        colorbar
        clim([1e-5 5e6]);
        colormap(colorMap_MP)
        colorbar( 'off' ) 
    end
end
