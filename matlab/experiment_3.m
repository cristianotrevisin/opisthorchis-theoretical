% In here we want to test the effect of fish market by running a simulation

clc; close all; clearvars;

% read OCN
OCN = build_OCN("OCN_A.mat",30*10000*10000);
par = common_parameters();
Time = 1:500*365;

set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

colorMap_MP = [ones(256,1)';linspace(1,0,256);linspace(1,0,256)]';



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

        setup = build_setup(OCN,par,33*800*1000,'seed',sd);
        setup.T = diag(min(setup.par.U*setup.H,setup.chi.*setup.F));
        setup.period = ones(length(Time),1);
    
        y0 = zeros(OCN.nNodes,4);
        y0(seeding_node,:) = [100/setup.H(seeding_node) 0 0 0];
    
        y = model_ODE(Time,setup.par,setup,y0');
    
        WH = y(:,1:4:end)';
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
for i = 1:3
    for ocnmap = 1:3
        counter = counter+1;
        if ocnmap == 1; WH_MAP = WH_1; elseif ocnmap == 2; WH_MAP = WH_2; else; WH_MAP = WH_3; end;
        subplot(3,3,counter)
        hold on
        loglog(Time/365,squeeze(WH_MAP(i,:,:)),'color',[0 0 0 0.1])
        loglog(Time/365,squeeze(WH_MAP(i,outlet_vec(ocnmap),:)),'color','#e4572e','LineWidth',1)
        loglog(Time/365,squeeze(WH_MAP(i,most_distant_node_vec(ocnmap),:)),'color','#f3a712','LineWidth',1)
        loglog(Time/365,squeeze(WH_MAP(i,seeding_node_vec(ocnmap),:)),'color','#29335c','LineWidth',1)
        set(gca,'XScale','log','YScale','log')
        xlim([1 Time(end)/365])
        ylim([1e-6 5e4])

        set(gca,'YTick',[1e-4 1e-2 1e0 1e2])
        set(gca,'XTick',[1e0 1e1 1e2])
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



%%
Times_Plot_Map = 365*[2 5 10 20];
OCN = build_OCN("OCN_A.mat",30*10000*10000);
counter = 0;
% figure
for i = 1:length(lambdas)
    for j = 1:length(Times_Plot_Map)
        counter = counter+1;
        %subplot(length(lambdas),length(Times_Plot_Map),counter)
        figure
        draw_OCN(OCN,squeeze(WH_1(i,:,Times_Plot_Map(j)))')
        if j == 1
            plot(OCN.geometry.SCX(seeding_node_vec(1))/OCN.cellsize,OCN.geometry.SCY(seeding_node_vec(1))/OCN.cellsize,'X','MarkerEdgeColor','red','MarkerSize',20)
        end
        set(gca,'ColorScale','log')
        colormap(colorMap_MP)
       colorbar
        clim([1e-6 5e4])
        set(gca,'FontSize',9)
    end
end

