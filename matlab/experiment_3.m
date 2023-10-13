% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;
% read OCN

OCN = build_OCN("OCN_A.mat",30*10000*10000);



par = common_parameters();

setup = build_setup(OCN,par,33*800*1000,'seed',3108);
setup.T = diag(1-setup.sigma);

outlet = find(OCN.SC_AccArea == max(OCN.SC_AccArea));
seeding_node = find(OCN.distW(outlet,:)==max(OCN.distW(outlet,:)));

y0 = zeros(OCN.nNodes,4);
y0(seeding_node,:) = [1/setup.H(seeding_node) 0 0 0];

Time = 1:50*365;

setup.period = ones(length(Time),1);


Times_Plot_Fig = [1 183 365 2*365 3*365 5*365 10*365 15*365 20*365 length(Time)];
Times_Plot_Map = 365*[2 3 5 10];
colorMap_IN = [linspace(016/256, 179/256,11); ...
    linspace(101/256, 021/256,11); ...
    linspace(171/256, 041/256,11)]';

colorMap_MP = [ones(256,1)';linspace(1,0,256);linspace(1,0,256)]';

% Only downstream
setup.par.lambda_FD = 0.01;
setup.par.lambda_FU = 0;

y = model_ODE(Time,setup.par,setup,y0');
WH_D = y(:,1:4:end);


% Both downstream and upstream
setup.par.lambda_FD = 0.01;
setup.par.lambda_FU = 0.01;

y = model_ODE(Time,setup.par,setup,y0');
WH_U = y(:,1:4:end);

%%

for tt = 1:length(Times_Plot_Map)
        figure
        draw_OCN(OCN,WH_D(Times_Plot_Map(tt),:)','Borders_Color','black')

        set(gca,'ColorScale','log')
        colorbar
        clim([1e-5 5e6]);
        colormap(colorMap_MP)
        colorbar( 'off' ) 
        if tt == 1
        plot(OCN.geometry.SCX(seeding_node)/OCN.cellsize,...
            OCN.geometry.SCY(seeding_node)/OCN.cellsize,...
            '.r','MarkerSize',10,'Marker','diamond')
        end
end


for tt = 1:length(Times_Plot_Map)
        figure
        draw_OCN(OCN,WH_U(Times_Plot_Map(tt),:)','Borders_Color','black')

        set(gca,'ColorScale','log')
        colorbar
        clim([1e-5 5e6]);
        colormap(colorMap_MP)

        if tt == 1
        plot(OCN.geometry.SCX(seeding_node)/OCN.cellsize,...
            OCN.geometry.SCY(seeding_node)/OCN.cellsize,...
            '.r','MarkerSize',10,'Marker','diamond')
        end

end

%%
for ocnmap = 1:3
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

    setup = build_setup(OCN,par,33*800*1000,'seed',sd);
    setup.T = diag(1-setup.sigma);
    setup.period = ones(length(Time),1);

    outlet = find(OCN.SC_AccArea == max(OCN.SC_AccArea));
    seeding_node = find(OCN.distW(outlet,:)==max(OCN.distW(outlet,:)));

    y0 = zeros(OCN.nNodes,4);
    y0(seeding_node,:) = [1/setup.H(seeding_node) 0 0 0];

    setup.par.lambda_FD = 0.01;
    setup.par.lambda_FU = 0;
    
    y_D = model_ODE(Time,setup.par,setup,y0');

    setup.par.lambda_FU = 0.01;
    
    y_U = model_ODE(Time,setup.par,setup,y0');

    if ocnmap == 1
        DIST1 = OCN.distW(seeding_node,:);

        WH1_D = y_D(:,1:4:end).*setup.H';
        WH1_U = y_U(:,1:4:end).*setup.H';
    elseif ocnmap == 2
        DIST2 = OCN.distW(seeding_node,:);

        WH2_D = y_D(:,1:4:end).*setup.H';
        WH2_U = y_U(:,1:4:end).*setup.H';
    elseif ocnmap == 3
        DIST3 = OCN.distW(seeding_node,:);

        WH3_D = y_D(:,1:4:end).*setup.H';
        WH3_U = y_U(:,1:4:end).*setup.H';
    end
end



% PLOT
% Mix distances
DIST = [DIST1 DIST2 DIST3];
WH_D = [DIST;
        [squeeze(WH1_D(Times_Plot_Fig,:))';
         squeeze(WH2_D(Times_Plot_Fig,:))';
         squeeze(WH3_D(Times_Plot_Fig,:))']'];

WH_D(:,abs(WH_D(end,:))< 1e-10) = [];

WH_D = sortrows(WH_D');

WH_U = [DIST;
        [squeeze(WH1_U(Times_Plot_Fig,:))';
         squeeze(WH2_U(Times_Plot_Fig,:))';
         squeeze(WH3_U(Times_Plot_Fig,:))']'];

WH_U(:,abs(WH_U(end,:))< 1e-10) = [];

WH_U = sortrows(WH_U');

figure()
hold on
for tt = fliplr(1:length(Times_Plot_Fig))
    area(WH_D(:,1)/1000,exp(smoothdata(log10(WH_D(:,tt+1)),'smoothingfactor',1)),'FaceColor',colorMap_IN(tt,:))
end
set(gca,'Yscale','log')
xlabel('Distance to first infected node [km]')
ylabel('I^H')
legend(num2str(fliplr(Times_Plot_Fig)'/365))
ylim([1e-8 1e5])
xlim([0 2500])

figure()
hold on
for tt = fliplr(1:length(Times_Plot_Fig))
    area(WH_U(:,1)/1000,exp(smoothdata(log10(WH_U(:,tt+1)),'smoothingfactor',1)),'FaceColor',colorMap_IN(tt,:))
end
set(gca,'Yscale','log')
xlabel('Distance to first infected node [km]')
ylabel('I^H')
ylim([1e-8 1e5])
xlim([0 2500])