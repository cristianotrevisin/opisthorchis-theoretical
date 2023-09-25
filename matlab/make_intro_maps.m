close all
clearvars
clc

OCN_A = build_OCN("OCN_A.mat");
OCN_B = build_OCN("OCN_B.mat");
OCN_C = build_OCN("OCN_C.mat");
% Add attributes
OCN_A.thrA = 55*10000*10000; 
OCN_B.thrA = 55*10000*10000; 
OCN_C.thrA = 55*10000*10000; 

par = common_parameters();
par.dF = 5;
par.dS = 50;
par.lambda_FU=0;
par.lambda_FD=0;
par.lambda_ED=0;
par.D = 100000;

% Create colormap
% Population humaine
cmap_H = [linspace(1, 0.9137254901960784, 256);...
    linspace(1, 0.6862745098039216, 256); ...
    linspace(1, 0.6392156862745098, 256)]';

cmap_F =  [linspace(1, 0.23137254901960785, 256);...
    linspace(1, 0.25098039215686274, 256); ...
    linspace(1, 0.35294117647058826, 256)]';

cmap_R =  [1 1 1; 0.5411764705882353 0.6039215686274509 0.3568627450980392];

% OCN_A

[setup_A1,par] = build_setup(OCN_A,par,33*800*1000,'seed',3108);
[setup_A2,~] = build_setup(OCN_A,par,33*800*1000,'seed',3108,'DownstreamAccumulation',true);

% Elevations
figure
title('Topographic map')
[cmap,zlimits] = ttcmap(OCN_A.FD.Z,'cmap','laos');
draw_OCN(OCN_A,OCN_A.FD.Z,'Borders_Color','k'); colormap(cmap);
cb=colorbar;
cb.Position = cb.Position + 1e-10;

figure
title('Human population (random)')
draw_OCN(OCN_A,setup_A1.H);
cb=colorbar;
cb.Position = cb.Position + 1e-10;
set(gca,'ColorScale','log')
colormap(cmap_H)

figure
title('Human population (downstream acc.)')
draw_OCN(OCN_A,setup_A2.H);
cb=colorbar;
cb.Position = cb.Position + 1e-10;
set(gca,'ColorScale','log')
colormap(cmap_H)


figure
title('Rice paddies')
draw_OCN(OCN_A,OCN_A.RicePaddy,'cmap',cmap_R,'Borders_Color','k')
colormap(cmap_R)

figure
title('Fish population')
colormap jet
draw_OCN(OCN_A,setup_A1.KF)
set(gca,'ColorScale','log')
cb=colorbar;
cb.Position = cb.Position + 1e-10;
colormap(cmap_F)

%%



figure
[cmap,zlimits] = ttcmap(OCN_B.FD.Z,'cmap','laos');
draw_OCN(OCN_B,OCN_B.FD.Z,'Borders_Color','k'); colormap(cmap); colorbar;

figure
[cmap,zlimits] = ttcmap(OCN_C.FD.Z,'cmap','laos');
draw_OCN(OCN_C,OCN_C.FD.Z,'Borders_Color','k'); colormap(cmap); colorbar;

figure
colormap sky
draw_OCN(OCN_A,setup.H)
set(gca,'ColorScale','log')
subplot(2,3,4)

colormap sky
draw_OCN(OCN_A,setup_A1.H)
set(gca,'ColorScale','log')
subplot(2,3,2)
colormap sky
draw_OCN(OCN_B,setup.H)
set(gca,'ColorScale','log')
subplot(2,3,5)
[setup,par] = build_setup(OCN_B,par,33*800*1000,'seed',3108,'DownstreamAccumulation',true);
colormap sky
draw_OCN(OCN_B,setup.H)
set(gca,'ColorScale','log')
subplot(2,3,3)
[setup,par] = build_setup(OCN_C,par,33*800*1000,'seed',3108);
colormap sky
draw_OCN(OCN_C,setup.H)
set(gca,'ColorScale','log')
subplot(2,3,6)
[setup,par] = build_setup(OCN_C,par,33*800*1000,'seed',3108,'DownstreamAccumulation',true);
colormap sky
draw_OCN(OCN_C,setup.H)
set(gca,'ColorScale','log')