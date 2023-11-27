close all
clearvars
clc

set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

OCN_A = build_OCN("OCN_A.mat",55*10000*10000);
OCN_B = build_OCN("OCN_B.mat",55*10000*10000);
OCN_C = build_OCN("OCN_C.mat",55*10000*10000);

par = common_parameters();


% Create colormap
% Population humaine
cmap_H = [linspace(1, 0.9137254901960784, 256);...
    linspace(1, 0.6862745098039216, 256); ...
    linspace(1, 0.6392156862745098, 256)]';

% Population de poissons
cmap_F =  [linspace(1, 0.23137254901960785, 256);...
    linspace(1, 0.25098039215686274, 256); ...
    linspace(1, 0.35294117647058826, 256)]';

% Présence de cultivations de riz
cmap_R =  [1 1 1; 0.5411764705882353 0.6039215686274509 0.3568627450980392];

% OCN_A
setup_A1 = build_setup(OCN_A,par,33*800*1000,'seed',3108);
setup_A2 = build_setup(OCN_A,par,33*800*1000,'seed',3108,'DownstreamAccumulation',true);

% OCN_B
setup_B1 = build_setup(OCN_B,par,33*800*1000,'seed',2507);
setup_B2 = build_setup(OCN_B,par,33*800*1000,'seed',2507,'DownstreamAccumulation',true);

% OCN_C
setup_C1 = build_setup(OCN_C,par,33*800*1000,'seed',0808);
setup_C2 = build_setup(OCN_C,par,33*800*1000,'seed',0808,'DownstreamAccumulation',true);

% FIND EXTREMES
maxZ = max([max(OCN_A.FD.Z),max(OCN_B.FD.Z),max(OCN_C.FD.Z)]);

maxH = max([max(setup_A1.H),max(setup_A2.H),...
    max(setup_B1.H),max(setup_B2.H),...
    max(setup_C1.H),max(setup_C2.H)]);
minH = min([min(setup_A1.H),min(setup_A2.H),...
    min(setup_B1.H),min(setup_B2.H),...
    min(setup_C1.H),min(setup_C2.H)]);

maxF = max([max(setup_A1.F),max(setup_A2.F),...
    max(setup_B1.F),max(setup_B2.F),...
    max(setup_C1.F),max(setup_C2.F)]);
minF = min([min(setup_A1.F(setup_A1.F>0)),min(setup_A2.F(setup_A2.F>0)),...
    min(setup_B1.F(setup_B1.F>0)),min(setup_B2.F(setup_B2.F>0)),...
    min(setup_C1.F(setup_C1.F>0)),min(setup_C2.F(setup_C2.F>0))]);

%% Plot 

[cmap_elevation,zlimits] = ttcmap(OCN_A.FD.Z,'cmap','laos'); 
figure
tiledlayout(5,3)

% elevations
ax11 = nexttile;
draw_OCN(OCN_A,OCN_A.FD.Z,'Borders_Color','k'); colormap(cmap_elevation); 
clim([0 maxZ]);

ax12 = nexttile;
draw_OCN(OCN_B,OCN_B.FD.Z,'Borders_Color','k'); colormap(cmap_elevation); 
clim([0 maxZ]);

ax13 = nexttile;
draw_OCN(OCN_C,OCN_C.FD.Z,'Borders_Color','k'); colormap(cmap_elevation); 
clim([0 maxZ]);

% human population (random)

ax21 = nexttile;
draw_OCN(OCN_A,setup_A1.H);
set(gca,'ColorScale','log')
colormap(cmap_H)
clim([minH maxH])

ax22 = nexttile;
draw_OCN(OCN_B,setup_B1.H);
set(gca,'ColorScale','log')
colormap(cmap_H)
clim([minH maxH])

ax23 = nexttile;
draw_OCN(OCN_C,setup_C1.H);
set(gca,'ColorScale','log')
colormap(cmap_H)
clim([minH maxH])

% human population (downstream acc)

ax31 = nexttile;
draw_OCN(OCN_A,setup_A2.H);
set(gca,'ColorScale','log')
colormap(cmap_H)
clim([minH maxH])

ax32 = nexttile;
draw_OCN(OCN_B,setup_B2.H);
set(gca,'ColorScale','log')
colormap(cmap_H)
clim([minH maxH])

ax33 = nexttile;
draw_OCN(OCN_C,setup_C2.H);
set(gca,'ColorScale','log')
colormap(cmap_H)
clim([minH maxH])

% rice paddies

ax41 = nexttile;
draw_OCN(OCN_A,OCN_A.RicePaddy,'cmap',cmap_R,'Borders_Color','k')
colormap(cmap_R)

ax42 = nexttile;
draw_OCN(OCN_B,OCN_B.RicePaddy,'cmap',cmap_R,'Borders_Color','k')
colormap(cmap_R)

ax43 = nexttile;
draw_OCN(OCN_C,OCN_C.RicePaddy,'cmap',cmap_R,'Borders_Color','k')
colormap(cmap_R)

% fish population

ax51 = nexttile;
draw_OCN(OCN_A,setup_A1.F)
set(gca,'ColorScale','log')
colormap(cmap_F)
clim([minF maxF])

ax52 = nexttile;
draw_OCN(OCN_B,setup_B1.F)
set(gca,'ColorScale','log')
colormap(cmap_F)
clim([minF maxF])

ax53 = nexttile;
draw_OCN(OCN_C,setup_C1.F)
set(gca,'ColorScale','log')
colormap(cmap_F)
clim([minF maxF])

colormap(ax11,cmap_elevation); colormap(ax12,cmap_elevation); colormap(ax13,cmap_elevation)
colormap(ax21,cmap_H); colormap(ax22,cmap_H); colormap(ax23,cmap_H)
colormap(ax31,cmap_H); colormap(ax32,cmap_H); colormap(ax33,cmap_H)
colormap(ax41,cmap_R); colormap(ax42,cmap_R); colormap(ax43,cmap_R)
colormap(ax51,cmap_F); colormap(ax52,cmap_F); colormap(ax53,cmap_F)
% plot for reference
figure
draw_OCN(OCN_A,setup_A1.H);
set(gca,'ColorScale','log')
colormap(cmap_H)
clim([minH maxH])
colorbar

figure
draw_OCN(OCN_A,OCN_A.FD.Z,'Borders_Color','k'); colormap(cmap_elevation); 
clim([0 maxZ]);
colorbar

figure
draw_OCN(OCN_A,setup_A1.F)
set(gca,'ColorScale','log')
colormap(cmap_F)
clim([minF maxF])
colorbar 



%%
% PLOT_fish_surplus (quantity)
colorMap_SP = [linspace(1, 0.45098039215686275,256);...
    linspace(1, 0.19607843137254902, 256);...
    linspace(1, 0.5098039215686274, 256)]';

purple = [0.45098039215686275,0.19607843137254902,0.5098039215686274];
colorMap_IN = [ones(1,256); linspace(1,0,256); linspace(1,0,256)]';

%% 
purple = "#6A0DAD";
close all
for i = 1:2
    if i == 1
        stp = setup_A1;
        ocn = OCN_A;
    elseif i == 2
        stp = setup_A2;
        ocn = OCN_A;
    end
    figure;
    draw_OCN(ocn,NaN)
    set(gca,'ColorScale','log')

    for nn = 1:ocn.nNodes
        for mm = 1:ocn.nNodes
            if mm~=nn && stp.T(nn,mm)>0
                l=line([ocn.geometry.SCX(nn)/ocn.cellsize ocn.geometry.SCX(mm)/ocn.cellsize],...
                    [ocn.geometry.SCY(nn)/ocn.cellsize ocn.geometry.SCY(mm)/ocn.cellsize],...
                    'linewidth',1);
                l.Color=[0,0,0,(stp.T(nn,mm)/max(stp.T-diag(diag(stp.T)),[],'all')).^0.1];
            end
        end
    end
    stp.sigma = stp.sigma/max(stp.sigma); stp.sigma(stp.sigma<=0)=NaN; 
    stp.delta = stp.delta/max(stp.delta); stp.delta(stp.delta<=0)=NaN; 
    scatter(ocn.geometry.SCX/ocn.cellsize,ocn.geometry.SCY/ocn.cellsize,500*stp.sigma.^0.5,'filled','MarkerFaceColor',purple)
    hold on
    scatter(ocn.geometry.SCX/ocn.cellsize,ocn.geometry.SCY/ocn.cellsize,500*stp.delta.^0.5,'MarkerEdgeColor',purple,'LineWidth',1)
end


