% EFFECT OF FISH TRADING

clc; close all; clearvars;

set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

par = common_parameters();
par.lambda_F = 0;
Time = 1:200*365;
Ds = [5e4 5e5 5e6];

% USER SHOULD SELECT WHICH OCNMAP TO USE
ocnmap = 1;

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


for scd = 1:length(Ds)
    par.D = Ds(scd);
    setup = build_setup(OCN,par,33*800*1000,'seed',sd);
    surplus = setup.sigma;
    
    setup_alt = setup;
    setup_alt.T = (min(setup_alt.par.U*setup_alt.H,setup_alt.chi.*setup_alt.F));
    
    SN = find(surplus==max(surplus));
    

    % Check if local endemic equilibrium is sustained
    tem = compute_EE(setup_alt);
    while  tem(SN,1)==0
        SN = find(surplus == max(surplus(surplus<surplus(SN))));
        tem = compute_EE(setup_alt);
    end

    SN1=SN;
    y0 = zeros(OCN.nNodes,4);
    y0(SN,1) = 0.1;
    y = model_ODE(Time,setup.par,setup,y0');
    lnt = size(y,1);
    if scd == length(Ds)
        H_C = setup.H;
    end
    WH(1:lnt,:,scd,1) = y(:,1:3:end);

    %DOWNSTREAM ACCUMULATION

    setup = build_setup(OCN,par,33*800*1000,'seed',sd,'DownstreamAccumulation',true);
    surplus = setup.sigma;
    setup_alt = setup;
    setup_alt.T = (min(setup_alt.par.U*setup_alt.H,setup_alt.chi.*setup_alt.F));

    SN = find(surplus==max(surplus));
    
    

    % Check if local endemic equilibrium is sustained
    tem = compute_EE(setup_alt);
    while  tem(SN,1)==0
        SN = find(surplus == max(surplus(surplus<surplus(SN))));
        tem = compute_EE(setup_alt);
    end

    SN2=SN;
    y0 = zeros(OCN.nNodes,4);
    y0(SN,1) = 0.1;

    y = model_ODE(Time,setup.par,setup,y0');
    lnt = size(y,1);

    WH(1:lnt,:,scd,2) = y(:,1:3:end);

    if scd == length(Ds)
        H_D = setup.H;
    end
end


%% Generate figures
figure
t = tiledlayout(3,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:3
    for j = 1:2
        if j == 1
            H = H_C;
        else
            H = H_D;
        end
        WHtot = squeeze(WH(:,:,i,j)); WHtot = sum(WHtot.*H',2,'omitmissing')/sum(H);
        nexttile
        hold on
        loglog(Time/365,WH(:,:,i,j)','color',[0 0 0 .15])
        if j == 1
            loglog(Time/365,WH(:,SN1,i,j)','color','red','LineWidth',1)
        elseif j == 2
            loglog(Time/365,WH(:,SN2,i,j)','color','red','LineWidth',1)
        end
        loglog(Time/365,WHtot,'color','k','LineWidth',1)
        ylim([1e-4  8e2])
        xlim([1 Time(end)/365])
        set(gca,'YScale','log')
        set(gca,'YTick',[1e-4 1e-2 1e0 1e2])
        %set(gca,'XTick',[1e0 1e1 1e2])
        if i < 3
            set(gca,'XTickLabel',[])
        end
        if i == 2 && j == 1
            ylabel ('I^H [worms/person]')
        end
        if i == 3
            xlabel ('Time [years]')
        end
        if j == 3 && i == 1
            ylabel ('D = 50 km')
        elseif j == 3 && i ==2
            ylabel ('D = 500 km')
        elseif j == 3 && i ==3
            ylabel ('D = 5000 km')
        end
        if i == 1 && j==1
            xlabel("Random allocation", "VerticalAlignment","top")
        elseif i == 1 && j==2
            xlabel("Downstream accumulation", "VerticalAlignment","top")
        end
    end
end
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [17 10]);
set(findall(gcf,'-property','FontSize'),'FontSize',9)



%%
Times_Plot_Map = 365*[1 5];

colorMap_MP = [ones(256,1)';linspace(1,0,256);linspace(1,0,256)]';

figure("Units","centimeters","PaperSize",[17*3,12*3])
t = tiledlayout(1,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(findall(gcf,'-property','FontSize'),'FontSize',9)
for d = 2
    for tt = 1:length(Times_Plot_Map)
        nexttile
        draw_OCN(OCN,WH(Times_Plot_Map(tt),:,d,1)')
        set(gca,'ColorScale','log')
        %colorbar
        clim([1e-6 5e2])
        colormap(colorMap_MP)
        %colorbar

    end
    for tt = 1:length(Times_Plot_Map)
        nexttile
        draw_OCN(OCN,WH(Times_Plot_Map(tt),:,d,2)')
        set(gca,'ColorScale','log')
        %colorbar
        clim([1e-6 5e2])
        colormap(colorMap_MP)
        %colorbar

    end
end
