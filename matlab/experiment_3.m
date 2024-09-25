%EFFECT OF FISH MOBILITY IN CONNECTED RIVER NETWORKS

clc; close all; clearvars;

par = common_parameters();
Time = 1:200*365;

set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

colorMap_MP = [ones(256,1)';linspace(1,0,256);linspace(1,0,256)]';

% USER SHOULD SELECT WHETHER TO USE DOWNSTREAM ACCUMULATION
downstream_accumulation=false;


lambdas = [1e-3 1e-2 1e-1];
outlet_vec = []; seeding_node_vec = []; most_distant_node_vec = [];
for ocnmap = 1:3
    if ocnmap == 1
        OCN = build_OCN("OCN_A.mat",30*10000*10000);
        sd = 3108;
        WH_1 = zeros(3,OCN.nNodes,length(Time));
    elseif ocnmap == 2
        OCN = build_OCN("OCN_B.mat",30*10000*10000);
        sd = 2507;
        WH_2 = zeros(3,OCN.nNodes,length(Time));
    elseif ocnmap == 3
        OCN = build_OCN("OCN_C.mat",30*10000*10000);
        sd = 0808;
        WH_3 = zeros(3,OCN.nNodes,length(Time));
    end

    
    outlet = find(OCN.SC_AccArea == max(OCN.SC_AccArea));
    seeding_node = find(OCN.distW(outlet,:)==max(OCN.distW(outlet,:)));
    while OCN.SC_RicePaddy_Area(seeding_node)==0
        seeding_node = find(OCN.distW(outlet,:) == max(OCN.distW(outlet,OCN.distW(outlet,:)<OCN.distW(outlet,seeding_node))));
    end
    most_distant_node = find(OCN.distW(seeding_node,:)==max(OCN.distW(seeding_node,:)));
    outlet_vec = [outlet_vec outlet]; seeding_node_vec = [seeding_node_vec seeding_node]; most_distant_node_vec = [most_distant_node_vec most_distant_node];
    

for i = 1:length(lambdas)
        par.lambda_F = lambdas(i);

        setup = build_setup(OCN,par,33*800*1000,'seed',sd,'DownstreamAccumulation',downstream_accumulation);
        setup.T = diag(min(setup.par.U*setup.H,setup.chi.*setup.F));
        setup.period = ones(length(Time),1);
    
        y0 = zeros(OCN.nNodes,4);
        y0(seeding_node,:) = [0.1 0 0 0];
    
        y = model_ODE(Time,setup.par,setup,y0');
    
        WH = y(:,1:3:end)';
        if ocnmap==1
            WH_1(i,:,:) = WH;
        elseif ocnmap == 2
            WH_2(i,:,:) = WH;
        else
            WH_3(i,:,:) = WH;
        end
    end
end
%%

counter = 0;
figure();
t = tiledlayout(3,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:3
    for ocnmap = 1:3
        counter = counter+1;
        if ocnmap == 1; WH_MAP = WH_1; elseif ocnmap == 2; WH_MAP = WH_2; else; WH_MAP = WH_3; end;
        nexttile()
        hold on
        loglog(Time/365,squeeze(WH_MAP(i,:,:)),'color',[0 0 0 0.1])
        loglog(Time/365,squeeze(WH_MAP(i,outlet_vec(ocnmap),:)),'color','#219ebc','LineWidth',1)
        loglog(Time/365,squeeze(WH_MAP(i,most_distant_node_vec(ocnmap),:)),'color','#ffb703','LineWidth',1)
        loglog(Time/365,squeeze(WH_MAP(i,seeding_node_vec(ocnmap),:)),'color','red','LineWidth',1)
        set(gca,'YScale','log')
        xlim([1 Time(end)/365])
        ylim([1e-6 2e2])

        set(gca,'YTick',[1e-4 1e-2 1e0 1e2])
        %set(gca,'XTick',[1e0 1e1 1e2])
        if ocnmap > 1
            set(gca,'YTickLabel',[])
        end
        if i < 3
            set(gca,'XTickLabel',[])
        end
      
        if ocnmap == 1 && i == 2
            ylabel ('I^H [worms/person]')
        elseif ocnmap == 2 && i == 3
            xlabel ('Time [years]')
        end
        if ocnmap == 3 && i == 1
            ylabel ('\lambda_F = 0.001 [1/day]')
        elseif ocnmap == 3 && i ==2
            ylabel ('\lambda_F = 0.01 [1/day]')
        elseif ocnmap == 3 && i ==3
            ylabel ('\lambda_F = 0.1 [1/day]')
        end
        set(gca,'FontSize',9)
    end
end
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [17 12]);
set(findall(gcf,'-property','FontSize'),'FontSize',9)


%% GENERATE FIGURES
Times_Plot_Map = 365*[1 5 10 50];
OCN = build_OCN("OCN_A.mat",30*10000*10000);
figure;
t = tiledlayout(3,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(findall(gcf,'-property','FontSize'),'FontSize',9)
for i = 1:length(lambdas)
    for j = 1:length(Times_Plot_Map)
        %subplot(length(lambdas),length(Times_Plot_Map),counter)
        nexttile()
        draw_OCN(OCN,squeeze(WH_1(i,:,Times_Plot_Map(j)))')
        set(gca,'ColorScale','log')
        colormap(colorMap_MP)
       %acolorbar
        clim([1e-6 5e2])
        set(gca,'FontSize',9)
    end
end

OCN = build_OCN("OCN_B.mat",30*10000*10000);
figure;
t = tiledlayout(3,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(findall(gcf,'-property','FontSize'),'FontSize',9)
for i = 1:length(lambdas)
    for j = 1:length(Times_Plot_Map)
        %subplot(length(lambdas),length(Times_Plot_Map),counter)
        nexttile()
        draw_OCN(OCN,squeeze(WH_2(i,:,Times_Plot_Map(j)))')
        set(gca,'ColorScale','log')
        colormap(colorMap_MP)
       %acolorbar
        clim([1e-6 5e2])
        set(gca,'FontSize',9)
    end
end

OCN = build_OCN("OCN_C.mat",30*10000*10000);
figure;
t = tiledlayout(3,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(findall(gcf,'-property','FontSize'),'FontSize',9)
for i = 1:length(lambdas)
    for j = 1:length(Times_Plot_Map)
        %subplot(length(lambdas),length(Times_Plot_Map),counter)
        nexttile()
        draw_OCN(OCN,squeeze(WH_3(i,:,Times_Plot_Map(j)))')
        set(gca,'ColorScale','log')
        colormap(colorMap_MP)
       %acolorbar
        clim([1e-6 5e2])
        set(gca,'FontSize',9)
    end
end