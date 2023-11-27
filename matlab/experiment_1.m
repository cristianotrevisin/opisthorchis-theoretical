% In here we want to investigate the existance of a local equilibrium


clc; close all; clearvars;

set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

par = common_parameters();

Time = 1:500*365; 

plt_RA = zeros(length(Time),4,3);
plt_DA = zeros(length(Time),4,3);
plt_UNI = zeros(length(Time),4,3);
p_eqi_RA = zeros(4,3);
p_eqi_DA = zeros(4,3);
p_eqi_UNI = zeros(4,3);

plot_maps = [1 0 0];
run_maps = 1%[1 2 3];
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

    setup_RA.T = diag(min(setup_RA.par.U*setup_RA.H,setup_RA.chi.*setup_RA.F));
    setup_DA.T = diag(min(setup_DA.par.U*setup_DA.H,setup_DA.chi.*setup_DA.F));


    % Now check endemic equilibrium in nodes (random allocation)
    eqi_RA = zeros(OCN.nNodes,4);
    for nn = 1:OCN.nNodes
        tmp = max(find_EE(setup_RA.par,setup_RA.H(nn), setup_RA.S(nn), setup_RA.A(nn),...
            setup_RA.chi(nn),setup_RA.par.epsilon(nn),setup_RA.par.theta(nn),...
            setup_RA.par.xi(nn),setup_RA.T(nn,nn)));
        if  isempty(tmp) || tmp(1)==0
            eqi_RA(nn,:)=0;
        else
            eqi_RA(nn,:)=tmp;
        end
    end
    EE_reached_RA = eqi_RA(:,1)>0;
    sum(EE_reached_RA)
    
    

    % Now check endemic equilibrium in nodes (downstream accumulation)
    eqi_DA = zeros(OCN.nNodes,4);
    for nn = 1:OCN.nNodes
        tmp = max(find_EE(setup_DA.par,setup_DA.H(nn), setup_DA.S(nn), setup_DA.A(nn),...
            setup_DA.chi(nn),setup_DA.par.epsilon(nn),setup_DA.par.theta(nn),...
            setup_DA.par.xi(nn),setup_DA.T(nn,nn)));
        if  isempty(tmp) || tmp(1)==0
            eqi_DA(nn,:)=0;
        else
            eqi_DA(nn,:)=tmp;
        end
    end
    EE_reached_DA = eqi_DA(:,1)>0;
    sum(EE_reached_DA)

    

    % Check endemic equilibrium in whole system
    eqi = find_EE(setup_UNI.par,setup_UNI.H,setup_UNI.S,setup_UNI.A,... 
        setup_UNI.chi,setup_UNI.par.epsilon,setup_UNI.par.theta,...
        setup_UNI.par.xi,setup_UNI.T);

    p_eqi_UNI(:,ocnmap) = max(eqi)';
    
    % Running simulations
    
    % Random allocation
    y0_RA = zeros(4,OCN.nNodes);
    y0_RA(1,:) = 100./setup_RA.H';
    
    % Downstream accumulation
    y0_DA = zeros(4,OCN.nNodes);
    y0_DA(1,:) = 100./setup_DA.H';
    
    % Unified
    y0_UNI = zeros(4,1);
    y0_UNI(1) = 100*setup_DA.nNodes./setup_UNI.H;


    % Launch simulations
    setup_RA.period = ones(length(Time),1); 
    setup_DA.period = ones(length(Time),1); setup_UNI.period = ones(length(Time),1);
    y_RA = model_ODE(Time,setup_RA.par,setup_RA,y0_RA);
    y_DA = model_ODE(Time,setup_DA.par,setup_DA,y0_DA);
    y_UNI = model_ODE(Time,setup_UNI.par,setup_UNI,y0_UNI);

    setup_RA.A(isnan(setup_RA.A)) = 0; setup_DA.A(isnan(setup_DA.A)) = 0;

    p_eqi_RA(1,ocnmap) = eqi_RA(:,1)'*setup_RA.H/sum(setup_RA.H); 
    p_eqi_RA(2,ocnmap) = eqi_RA(:,2)'*setup_RA.A/sum(setup_RA.A); 
    p_eqi_RA(3,ocnmap) = eqi_RA(:,3)'*setup_RA.S/sum(setup_RA.S);
    p_eqi_RA(4,ocnmap) = eqi_RA(:,4)'*setup_RA.F/sum(setup_RA.F);


    p_eqi_DA(1,ocnmap) = eqi_DA(:,1)'*setup_DA.H/sum(setup_DA.H); 
    p_eqi_DA(2,ocnmap) = eqi_DA(:,2)'*setup_DA.A/sum(setup_DA.A); 
    p_eqi_DA(3,ocnmap) = eqi_DA(:,3)'*setup_DA.S/sum(setup_DA.S); 
    p_eqi_DA(4,ocnmap) = eqi_DA(:,4)'*setup_DA.F/sum(setup_DA.F);


    % Agglomerate results
    plt_RA(:,1,ocnmap) = y_RA(:,1:4:end)*setup_RA.H/(sum(setup_RA.H));
    plt_RA(:,2,ocnmap) = y_RA(:,2:4:end)*setup_RA.A/(sum(setup_RA.A,'omitmissing'));
    plt_RA(:,3,ocnmap) = y_RA(:,3:4:end)*setup_RA.S/(sum(setup_RA.S));
    plt_RA(:,4,ocnmap) = y_RA(:,4:4:end)*setup_RA.F/(sum(setup_RA.F));
    
    WH_RA = y_RA(:,1:4:end);
    WH_DA = y_DA(:,1:4:end);
    
    plt_DA(:,1,ocnmap) = y_DA(:,1:4:end)*setup_DA.H/(sum(setup_DA.H));
    plt_DA(:,2,ocnmap) = y_DA(:,2:4:end)*setup_DA.A/(sum(setup_DA.A,'omitmissing'));
    plt_DA(:,3,ocnmap) = y_DA(:,3:4:end)*setup_DA.S/(sum(setup_DA.S));
    plt_DA(:,4,ocnmap) = y_DA(:,4:4:end)*setup_DA.F/(sum(setup_DA.F));

    plt_UNI(:,1,ocnmap) = y_UNI(:,1);
    plt_UNI(:,2,ocnmap) = y_UNI(:,2);
    plt_UNI(:,3,ocnmap) = y_UNI(:,3);
    plt_UNI(:,4,ocnmap) = y_UNI(:,4);

    if plot_maps(ocnmap) == 1

        figure
        if any(EE_reached_RA==1)
            loglog(Time/365,WH_RA(:,EE_reached_RA==1),'color','#cc2936');
        end
        hold on
        if any(EE_reached_RA==0)
            loglog(Time/365,WH_RA(:,EE_reached_RA==0),'color','#08415c');
        end
        xlim([1/10 Time(end)/365])
        ylim([1/max(setup_RA.H) max(WH_RA,[],'all')])
        ylabel('I^H')
        xlabel('Years')


        figure
        draw_OCN(OCN,EE_reached_RA,'binary',true)


        figure
        if any(EE_reached_DA==1)
            loglog(Time/365,WH_DA(:,EE_reached_DA==1),'color','#cc2936');
        end
        hold on
        if any(EE_reached_DA==0)
            loglog(Time/365,WH_DA(:,EE_reached_DA==0),'color','#08415c');
        end
        xlim([1/10 Time(end)/365])
        ylim([1/max(setup_DA.H) max(WH_DA,[],'all')])
        ylabel('I^H')
        xlabel('Years')
        set(gca,'yscale','log','xscale','log')
    
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

for ocnmap = 1:3
        ind_cut = min([find(plt_RA(:,1,ocnmap)<0,1,'first'),...
            find(plt_RA(:,2,ocnmap)<0,1,'first'),...
            find(plt_RA(:,3,ocnmap)<0,1,'first'),...
            find(plt_RA(:,4,ocnmap)<0,1,'first')]);
        plt_RA(ind_cut:end,:,ocnmap) = NaN;
        ind_cut = min([find(plt_DA(:,1,ocnmap)<0,1,'first'),...
            find(plt_DA(:,2,ocnmap)<0,1,'first'),...
            find(plt_DA(:,3,ocnmap)<0,1,'first'),...
            find(plt_DA(:,4,ocnmap)<0,1,'first')]);
        plt_DA(ind_cut:end,:,ocnmap) = NaN;
end


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
    plot(Time/365,ones(length(Time),1)*max(p_eqi_RA(1,ocnmap)),'--','color','#ef476f')
    plot(Time/365,ones(length(Time),1)*max(p_eqi_DA(1,ocnmap)),'--','color','#ffd166')
end
ylabel('I^H [worms/person]')
set(gca,'XTick',[0.01 0.1 1 10 100])
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
    plot(Time/365,ones(length(Time),1)*max(p_eqi_RA(2,ocnmap)),'--','color','#ef476f')
    plot(Time/365,ones(length(Time),1)*max(p_eqi_DA(2,ocnmap)),'--','color','#ffd166')
end
ylabel('E [eggs/m^3]')
set(gca,'XTick',[0.01 0.1 1 10 100])
xlim([Time(1)/365 Time(end)/365])
set(gca,'XTickLabel',[])
set(gca,'XScale','log','YScale','log')
box off
grid minor

nexttile()
hold on
for ocnmap = 1:3
    plot(Time/365,plt_UNI(:,3,ocnmap),'color','#26547c')
    plot(Time/365,plt_RA(:,3,ocnmap),'color','#ef476f')
    plot(Time/365,plt_DA(:,3,ocnmap),'color','#ffd166')
    plot(Time/365,ones(length(Time),1)*max(p_eqi_UNI(3,ocnmap)),'--','color','#26547c')
    plot(Time/365,ones(length(Time),1)*max(p_eqi_RA(3,ocnmap)),'--','color','#ef476f')
    plot(Time/365,ones(length(Time),1)*max(p_eqi_DA(3,ocnmap)),'--','color','#ffd166')
end
ylabel('Y^S')
set(gca,'XTick',[0.01 0.1 1 10 100])
xlim([Time(1)/365 Time(end)/365])
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
    plot(Time/365,ones(length(Time),1)*max(p_eqi_RA(4,ocnmap)),'--','color','#ef476f')
    plot(Time/365,ones(length(Time),1)*max(p_eqi_DA(4,ocnmap)),'--','color','#ffd166')
end
ylabel('I^F [cysts/fish]')
set(gca,'XTick',[0.01 0.1 1 10 100])
xlim([Time(1)/365 Time(end)/365])
set(gca,'XScale','log','YScale','log')
xlabel('Years')
box off
grid minor