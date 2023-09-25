% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;
% read OCN

OCN = build_OCN("OCN_A.mat");
% Add attributes
OCN.thrA = 30*10000*10000; 


par = common_parameters();
par.dF = 10;
par.dS = 30;
par.c = 3.0280e-08;

par.D = 100000;

[setup] = build_setup(OCN,par,33*800*1000,'seed',3108);
% Add parameter after volume
betaHS = 9.16e-11;
par.beta_E = betaHS*mean(setup.V)*par.mu_E/(par.rho_E-mean(setup.V)*betaHS*3102/7.5/0.5);



% PLOT_fish_deficit (quantity)
par.U = 1;
Surplus = par.c*setup.H.*setup.KF - par.U*setup.H; Surplus(Surplus<0)=0;
LC = par.c*setup.H.*setup.KF; LC = repmat(LC',OCN.nNodes,1);
TRA = setup.T.*LC;
figure
colorMap_SP = [linspace(1, 0.45098039215686275,256);...
    linspace(1, 0.19607843137254902, 256);...
    linspace(1, 0.5098039215686274, 256)]';
colormap(colorMap_SP);

colorMap_IN = [ones(1,256); linspace(1,0,256); linspace(1,0,256)]';
draw_OCN(OCN,Surplus,'Borders_Color','black')
%draw_OCN(OCN,Surplus==0,'Borders_Color','black')
set(gca,'ColorScale','log')
colorbar 
% for nn = 1:OCN.nNodes
%     for mm = 1:OCN.nNodes
%         if mm~=nn && TRA(nn,mm)>0
%             l=line([OCN.geometry.SCX(nn)/OCN.cellsize OCN.geometry.SCX(mm)/OCN.cellsize],...
%                 [OCN.geometry.SCY(nn)/OCN.cellsize OCN.geometry.SCY(mm)/OCN.cellsize],...
%                 'linewidth',1);%2*(TRA(nn,mm)/max(TRA-diag(diag(TRA)),[],'all'))^0.25);
%             l.Color=[0,0,0,(TRA(nn,mm)/max(TRA-diag(diag(TRA)),[],'all')).^0.4];
%         end
%     end
% end
% for sc = 1:OCN.nNodes
%     plot(OCN.geometry.SCX(sc)/OCN.cellsize,OCN.geometry.SCY(sc)/OCN.cellsize,'.r','MarkerSize',0.5+1.5*log(setup.H(sc)))
%     %text(OCN.geometry.SCX(sc)/OCN.cellsize,OCN.geometry.SCY(sc)/OCN.cellsize,num2str(sc),'Color','k')
% end
%%
surplus = setup.sigma.*setup.KF;

y0 = zeros(OCN.nNodes,5);
y0(:,4) = setup.KF;
y0(Surplus==max(Surplus),:) = [1 0 0 setup.KF(Surplus==max(Surplus)) 0];

Time = 1:10*365;
setup.nNodes = OCN.nNodes;

par.lambda_FU=0;
par.lambda_FD=0;
par.lambda_ED=0;

y = model_ODE(Time,par,setup,y0');

WH = y(:,1:5:end);
figure()
semilogy(Time/365,WH);
for nn = 1:OCN.nNodes
    text(Time(end)/365,WH(end,nn),num2str(nn));
end

%%
close all
figure()
draw_OCN(OCN,WH(end,:)','Borders_Color','black')
%draw_OCN(OCN,Surplus==0,'Borders_Color','black')
set(gca,'ColorScale','log')
colorbar
clim([1e-5 5e6]);
%colorbar( 'off' ) 
colormap(colorMap_IN)

%% TEST DIFFERENT D
Ds = [1e3 1e4 1e5 1e6];
WHs = zeros(length(Ds),length(Time));
Times = [1 30 365 10*365];
figure
cntr = 0;
for i = 1:length(Ds)
    par.D = Ds(i);
    [setup] = build_setup(OCN,par,33*800*1000,'seed',3108);
    surplus = setup.sigma.*setup.KF;
    y0 = zeros(OCN.nNodes,5);
    y0(:,4) = setup.KF;
    y0(Surplus==max(Surplus),:) = [1 0 0 setup.KF(Surplus==max(Surplus)) 0];
    
    
    y = model_ODE(Time,par,setup,y0');
    WH = y(:,1:5:end);
    WHs(i,:) = WH*setup.H/sum(setup.H);

    for tt = 1:length(Times)
        cntr = cntr+1;
        subplot(4,4,cntr)
        draw_OCN(OCN,WH(Times(tt),:)','Borders_Color','black')
        set(gca,'ColorScale','log')
        colorbar
        clim([1e-5 5e6]);
        colormap(colorMap_IN)
        colorbar( 'off' ) 
    end



end

figure
semilogy(Time/365,WHs)

