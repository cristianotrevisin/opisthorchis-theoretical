% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;

set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

par = common_parameters();
Time = 1:500*365;
Ds = [1e4 1e5 1e6];

OCN = build_OCN("OCN_A.mat",30*10000*10000);
sd = 3108;




for scd = 1:length(Ds)
    par.D = Ds(scd);
    setup = build_setup(OCN,par,33*800*1000,'seed',sd);
    setup.period = ones(length(Time),1);
    surplus = setup.sigma;

    y0 = zeros(OCN.nNodes,4);
    y0(surplus==max(surplus),1) = 100./setup.H(surplus==max(surplus));
    
    y = model_ODE(Time,setup.par,setup,y0');
    if scd == length(Ds)
        H_C = setup.H;
    end
    WH(:,:,scd,1) = y(:,1:4:end);

    setup = build_setup(OCN,par,33*800*1000,'seed',sd,'DownstreamAccumulation',true);
    setup.period = ones(length(Time),1);
    surplus = setup.sigma;

    y0 = zeros(OCN.nNodes,4);
    y0(surplus==max(surplus),1) = 100./setup.H(surplus==max(surplus));

    y = model_ODE(Time,setup.par,setup,y0');

    WH(:,:,scd,2) = y(:,1:4:end);

    if scd == length(Ds)
        H_D = setup.H;
    end
end


%%
figure
for i = 1:3
    for j = 1:2
        if j == 1
            H = H_C;
        else
            H = H_D;
        end
        WHtot = squeeze(WH(:,:,i,j)); WHtot = sum(WHtot.*H',2,'omitmissing')/sum(H);
        if j == 1
            subplot(3,2,2*i-1)
        else
            subplot(3,2,2*i)
        end
        hold on
        loglog(Time/365,WH(:,:,i,j)','color',[0 0 0 .15])
        loglog(Time/365,WH(:,surplus==max(surplus),i,j)','color','red','LineWidth',1)
        loglog(Time/365,WHtot,'color','k','LineWidth',1)
        ylim([1e-8 1e4])
        xlim([1 Time(end)/365])
        set(gca,'XScale','log','YScale','log')
        set(gca,'YTick',[1e-6 1e-3 1e0 1e3])
        set(gca,'XTick',[1e0 1e1 1e2])
        if j > 1
            set(gca,'YTickLabel',[])
        end
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
            ylabel ('D = 10 km')
        elseif j == 3 && i ==2
            ylabel ('D = 100 km')
        elseif j == 3 && i ==3
            ylabel ('D = 1000 km')
        end
    end
end



%%
Times_Plot_Map = 365*[5 10 50 100];


OCN = build_OCN("OCN_A.mat",30*10000*10000);
        sd = 3108;
colorMap_MP = [ones(256,1)';linspace(1,0,256);linspace(1,0,256)]';

cntr = 0;
for d = 2
    for tt = 1:length(Times_Plot_Map)
        figure;
        draw_OCN(OCN,WH(Times_Plot_Map(tt),:,d,2)')
        set(gca,'ColorScale','log')
        colorbar
        clim([1e-8 1e4])
        colormap(colorMap_MP)
        colorbar('off')

    end
end
