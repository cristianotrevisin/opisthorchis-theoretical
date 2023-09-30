% In here we want to investigate the existance of a local equilibrium


clc; close all; clearvars;

par = common_parameters();

Time = 1:1000*365; 

plt_RA = zeros(length(Time),4,3);
plt_DA = zeros(length(Time),4,3);
plt_UNI = zeros(length(Time),4,3);
p_eqi_RA = zeros(4,3);
p_eqi_DA = zeros(4,3);
p_eqi_UNI = zeros(4,3);

plot_maps = [0 0 0];
run_maps = [1 2 3];
for ocnmap = run_maps
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

    % Load setups

    setup_RA = build_setup(OCN,par,33*800*1000,'seed',sd);
    setup_DA = build_setup(OCN,par,33*800*1000,'seed',sd,'DownstreamAccumulation',true);
    setup_UNI = build_setup(OCN,par,33*800*1000,'seed',sd,'Unify',true);


    % Prevent transportation of raw fish outside community
    setup_RA.T = eye(OCN.nNodes); %setup.T(eye(OCN.nNodes)==1) = 1-setup.sigma;
    setup_DA.T = eye(OCN.nNodes); %setup.T(eye(OCN.nNodes)==1) = 1-setup.sigma;



    % Now check endemic equilibrium in nodes (random allocation)
    eqi_RA = zeros(OCN.nNodes,4);
    for nn = 1:OCN.nNodes
        tmp = max(find_EE(setup_RA.par,setup_RA.par.c*setup_RA.H(nn),setup_RA.H(nn),...
            setup_RA.S(nn),setup_RA.F(nn),setup_RA.V(nn)));
        if  isempty(tmp) || tmp(1)==0
            eqi_RA(nn,:)=0;
        else
            eqi_RA(nn,:)=tmp;
        end
    end
    EE_reached_RA = eqi_RA(:,1)>0;
    sum(EE_reached_RA)
    
    p_eqi_RA(1,ocnmap) = eqi_RA(:,1)'*setup_RA.H/sum(setup_RA.H); 
    p_eqi_RA(2,ocnmap) = eqi_RA(:,2)'*setup_RA.V/sum(setup_RA.V); 
    p_eqi_RA(3,ocnmap) = eqi_RA(:,3)'*setup_RA.S/sum(setup_RA.S);
    p_eqi_RA(4,ocnmap) = eqi_RA(:,4)'*setup_RA.F/sum(setup_RA.F);

    % Now check endemic equilibrium in nodes (downstream accumulation)
    eqi_DA = zeros(OCN.nNodes,4);
    for nn = 1:OCN.nNodes
        tmp = max(find_EE(setup_DA.par,setup_DA.par.c*setup_DA.H(nn),setup_DA.H(nn),...
            setup_DA.S(nn),setup_DA.F(nn),setup_DA.V(nn)));
        if  isempty(tmp) || tmp(1)==0
            eqi_DA(nn,:)=0;
        else
            eqi_DA(nn,:)=tmp;
        end
    end
    EE_reached_DA = ~isnan(eqi_DA(:,1));
    sum(EE_reached_DA)

    p_eqi_DA(1,ocnmap) = eqi_DA(:,1)'*setup_DA.H/sum(setup_DA.H); 
    p_eqi_DA(2,ocnmap) = eqi_DA(:,2)'*setup_DA.V/sum(setup_DA.V); 
    p_eqi_DA(3,ocnmap) = eqi_DA(:,3)'*setup_DA.S/sum(setup_DA.S); 
    p_eqi_DA(4,ocnmap) = eqi_DA(:,4)'*setup_DA.F/sum(setup_DA.F);

    % Check endemic equilibrium in whole system
    eqi = find_EE(setup_UNI.par,...
        setup_UNI.par.c*setup_UNI.H,...
        setup_UNI.H,setup_UNI.S,setup_UNI.F,setup_UNI.V);

    p_eqi_UNI(:,ocnmap) = max(eqi)';
    
    % Running simulations
    
    % Random allocation
    y0_RA = zeros(4,OCN.nNodes);
    y0_RA(1,:) = 1./setup_RA.H';
    
    % Downstream accumulation
    y0_DA = zeros(4,OCN.nNodes);
    y0_DA(1,:) = 1./setup_DA.H';
    
    % Unified
    y0_UNI = zeros(4,1);
    y0_UNI(1) = setup_DA.nNodes./setup_UNI.H;


    % Launch simulations
    setup_RA.period = ones(length(Time),1); 
    setup_DA.period = ones(length(Time),1); setup_UNI.period = ones(length(Time),1);
    y_RA = model_ODE(Time,setup_RA.par,setup_RA,y0_RA);
    y_DA = model_ODE(Time,setup_DA.par,setup_DA,y0_DA);
    y_UNI = model_ODE(Time,setup_UNI.par,setup_UNI,y0_UNI);


    % Agglomerate results
    plt_RA(:,1,ocnmap) = y_RA(:,1:4:end)*setup_RA.H/(sum(setup_RA.H));
    plt_RA(:,2,ocnmap) = y_RA(:,2:4:end)*setup_RA.V/(sum(setup_RA.V));
    plt_RA(:,3,ocnmap) = y_RA(:,3:4:end)*setup_RA.S/(sum(setup_RA.S));
    plt_RA(:,4,ocnmap) = y_RA(:,4:4:end)*setup_RA.F/(sum(setup_RA.F));
    
    WH_RA = y_RA(:,1:4:end);
    WH_DA = y_DA(:,1:4:end);
    
    plt_DA(:,1,ocnmap) = y_DA(:,1:4:end)*setup_DA.H/(sum(setup_DA.H));
    plt_DA(:,2,ocnmap) = y_DA(:,2:4:end)*setup_DA.V/(sum(setup_DA.V));
    plt_DA(:,3,ocnmap) = y_DA(:,3:4:end)*setup_DA.S/(sum(setup_DA.S));
    plt_DA(:,4,ocnmap) = y_DA(:,4:4:end)*setup_DA.F/(sum(setup_DA.F));

    plt_UNI(:,1,ocnmap) = y_UNI(:,1);
    plt_UNI(:,2,ocnmap) = y_UNI(:,2);
    plt_UNI(:,3,ocnmap) = y_UNI(:,3);
    plt_UNI(:,4,ocnmap) = y_UNI(:,4);

    if plot_maps(ocnmap) == 1

        figure
        loglog(Time/365,WH_RA(:,EE_reached_RA==1),'color','#cc2936');
        hold on
        loglog(Time/365,WH_RA(:,EE_reached_RA==0),'color','#08415c');
        xlim([1/10 Time(end)/365])
        ylim([1/max(setup_RA.H) max(WH_RA,[],'all')])
        ylabel('I^H')
        xlabel('Years')


        figure
        draw_OCN(OCN,EE_reached_RA,'binary',true)


        figure
        loglog(Time/365,WH_DA(:,EE_reached_DA==1),'color','#cc2936');
        hold on
        loglog(Time/365,WH_DA(:,EE_reached_DA==0),'color','#08415c');
        xlim([1/10 Time(end)/365])
        ylim([1/max(setup_DA.H) max(WH_DA,[],'all')])
        ylabel('I^H')
        xlabel('Years')
    
        figure
        draw_OCN(OCN,EE_reached_DA,'binary',true)
    end
end
%%

plt_UNI_mean = squeeze(mean(plt_UNI,3));
plt_RA_mean = squeeze(mean(plt_RA,3));
plt_DA_mean = squeeze(mean(plt_DA,3));

p_eqi_DA_mean = squeeze(mean(p_eqi_DA,2));
p_eqi_RA_mean = squeeze(mean(p_eqi_RA,2));
p_eqi_UNI_mean = squeeze(mean(p_eqi_UNI,2));
%%

figure()

tiledlayout(2,2)

nexttile()
hold on
plot(Time/365,plt_UNI_mean(:,1),'color','#26547c')
plot(Time/365,plt_RA_mean(:,1),'color','#ef476f')
plot(Time/365,plt_DA_mean(:,1),'color','#ffd166')
plot(Time/365,ones(length(Time),1)*max(p_eqi_UNI_mean(1)),'--','color','#26547c')
plot(Time/365,ones(length(Time),1)*max(p_eqi_RA_mean(1)),'--','color','#ef476f')
plot(Time/365,ones(length(Time),1)*max(p_eqi_DA_mean(1)),'--','color','#ffd166')
ylabel('I^H [w/h]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
set(gca,'XScale','log','YScale','log')
box off
grid minor 

nexttile()
hold on
plot(Time/365,plt_UNI_mean(:,2),'color','#26547c')
plot(Time/365,plt_RA_mean(:,2),'color','#ef476f')
plot(Time/365,plt_DA_mean(:,2),'color','#ffd166')
plot(Time/365,ones(length(Time),1)*max(p_eqi_UNI_mean(2)),'--','color','#26547c')
plot(Time/365,ones(length(Time),1)*max(p_eqi_RA_mean(2)),'--','color','#ef476f')
plot(Time/365,ones(length(Time),1)*max(p_eqi_DA_mean(2)),'--','color','#ffd166')
ylabel('E [e/m^3]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
set(gca,'XScale','log','YScale','log')
box off
grid minor

nexttile()
hold on
plot(Time/365,plt_UNI_mean(:,3),'color','#26547c')
plot(Time/365,plt_RA_mean(:,3),'color','#ef476f')
plot(Time/365,plt_DA_mean(:,3),'color','#ffd166')
plot(Time/365,ones(length(Time),1)*max(p_eqi_UNI_mean(3)),'--','color','#26547c')
plot(Time/365,ones(length(Time),1)*max(p_eqi_RA_mean(3)),'--','color','#ef476f')
plot(Time/365,ones(length(Time),1)*max(p_eqi_DA_mean(3)),'--','color','#ffd166')
ylabel('Y^S [%]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
set(gca,'XScale','log','YScale','log')
box off
grid minor 


nexttile()
hold on
plot(Time/365,plt_UNI_mean(:,4),'color','#26547c')
plot(Time/365,plt_RA_mean(:,4),'color','#ef476f')
plot(Time/365,plt_DA_mean(:,4),'color','#ffd166')
plot(Time/365,ones(length(Time),1)*max(p_eqi_UNI_mean(4)),'--','color','#26547c')
plot(Time/365,ones(length(Time),1)*max(p_eqi_RA_mean(4)),'--','color','#ef476f')
plot(Time/365,ones(length(Time),1)*max(p_eqi_DA_mean(4)),'--','color','#ffd166')
ylabel('I^F [m/f]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XScale','log','YScale','log')
xlabel('Years')
box off
grid minor
%%
figure()

tiledlayout(2,2)

nexttile()
hold on
for ocnmap = 1:3
    plot(Time/365,plt_UNI(:,1,ocnmap),'color','#26547c')
    plot(Time/365,plt_RA(:,1,ocnmap),'color','#ef476f')
    plot(Time/365,plt_DA(:,1,ocnmap),'color','#ffd166')
    plot(Time/365,ones(length(Time),1)*max(p_eqi_UNI(1,ocnmap)),'--','color','#26547c')
end
ylabel('I^H [w/h]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
set(gca,'XScale','log','YScale','log')
box off
grid minor 

nexttile()
hold on
for ocnmap = 1:3
    plot(Time/365,plt_UNI(:,2,ocnmap),'color','#26547c')
    plot(Time/365,plt_RA(:,2,ocnmap),'color','#ef476f')
    plot(Time/365,plt_DA(:,2,ocnmap),'color','#ffd166')
    plot(Time/365,ones(length(Time),1)*max(p_eqi_UNI(2,ocnmap)),'--','color','#26547c')
end
ylabel('E [e/m^3]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
set(gca,'XScale','log','YScale','log')
box off

nexttile()
hold on
for ocnmap = 1:3
    plot(Time/365,plt_UNI(:,3,ocnmap),'color','#26547c')
    plot(Time/365,plt_RA(:,3,ocnmap),'color','#ef476f')
    plot(Time/365,plt_DA(:,3,ocnmap),'color','#ffd166')
    plot(Time/365,ones(length(Time),1)*max(p_eqi_UNI(3,ocnmap)),'--','color','#26547c')
end
ylabel('Y^S [%]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
set(gca,'XScale','log','YScale','log')
box off
grid minor 
xlabel('Years')

nexttile()
hold on
for ocnmap = 1:3
    plot(Time/365,plt_UNI(:,4,ocnmap),'color','#26547c')
    plot(Time/365,plt_RA(:,4,ocnmap),'color','#ef476f')
    plot(Time/365,plt_DA(:,4,ocnmap),'color','#ffd166')
    plot(Time/365,ones(length(Time),1)*max(p_eqi_UNI(4,ocnmap)),'--','color','#26547c')
end
ylabel('I^F [m/f]')
xlim([Time(1)/365 Time(end)/365])
set(gca,'XScale','log','YScale','log')
xlabel('Years')
box off
grid minor