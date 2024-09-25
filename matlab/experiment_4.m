clc; close all; clearvars;
% read OCN
set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')



par = common_parameters();
Dvec = 5*10.^[2:0.25:6];
Lvec = 10.^[-4:0.25:-1];

ocnmap = 3;
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



for dv = 1:length(Dvec)
    dv
    par.D = Dvec(dv);
    for lv = 1:length(Lvec)
        lv
        par.lambda_F = Lvec(lv);

        setup = build_setup(OCN,par,33*800*1000,'seed',sd);
        outlet = find(OCN.SC_AccArea == max(OCN.SC_AccArea));
        SN = find(OCN.distW(outlet,:)==max(OCN.distW(outlet,:)));

        out = get_simulation_equilibrium(setup,SN);

        EQ1(dv,lv) = out(end);


        TM1(dv,lv) = find(out>0.95*out(end),1,'first')/365;
        out1(:,dv,lv) = out;


        setup = build_setup(OCN,par,33*800*1000,'seed',sd,'DownstreamAccumulation',true);
        outlet = find(OCN.SC_AccArea == max(OCN.SC_AccArea));
        SN = find(OCN.distW(outlet,:)==max(OCN.distW(outlet,:)));

        out = get_simulation_equilibrium(setup,SN);

        EQ2(dv,lv) = out(end);

        

        TM2(dv,lv) = find(out>0.95*out(end),1,'first')/365;

        out2(:,dv,lv) = out;
    end
end


%%
MINT = min(min(TM1(:)),min(TM2(:)));
MAXT = max(max(TM1(:)),max(TM2(:)));

MINE = min(min(EQ1(:)),min(EQ2(:)));
MAXE = max(max(EQ1(:)),max(EQ2(:)));



[X, Y] = meshgrid(Dvec/1000, Lvec);
figure
t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile()
contourf(X, Y, EQ1');
pcolor(X, Y, EQ1');
set(gca,'XScale','log','YScale','log')
shading interp
clim([MINE MAXE])
set(gca,'XTickLabel',[])
ylabel('\lambda_F [1/day]')


nexttile()
contourf(X, Y, EQ2');
pcolor(X, Y, EQ2');
set(gca,'XScale','log','YScale','log')
shading interp
clim([MINE MAXE]);
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
colorbar

nexttile()
contourf(X, Y, TM1');
pcolor(X, Y, TM1');
set(gca,'XScale','log','YScale','log')
shading interp
clim([MINT MAXT])
ylabel('\lambda_F [1/day]')
set(gca,'XTick',[1 10 100 1000 10000])
xlabel('D [km]')

nexttile()
contourf(X, Y, TM2');
pcolor(X, Y, TM2');
set(gca,'XScale','log','YScale','log')
shading interp
clim([MINT MAXT])
colorbar
set(gca,'YTickLabel',[])
xlabel('D [km]')
set(gca,'XTick',[1 10 100 1000 10000])

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [17 12]);
set(findall(gcf,'-property','FontSize'),'FontSize',9)
%%
function out = get_simulation_equilibrium(setup,SN)
    Time = 1:500*365;
    
    y0 = zeros(setup.nNodes,4);

    y0(SN,1) = 0.1;

    
    y = model_ODE_stiff(Time,setup.par,setup,y0');

    out = y(:,1:3:end)*setup.H/(sum(setup.H));

end


function y = model_ODE_stiff(Time,par,setup,y0)

    %node_out = find(sum(W,1)==0);

    y=odemodel(par,...
        setup.nNodes,... #number of nodes
        setup.H,... #human population
        setup.S,... #snail population
        setup.F,...#carrying fish population
        setup.T,... #caught fish trade
        setup.A,... #local area for snails
        setup.W,... #hydrological connectivity
        setup.chi,...#fish catch rate
        setup.par.xi,... #urbanization reduction of snails' exposure
        setup.par.epsilon,... #fraction of fish consumed raw
        setup.par.theta,... #fish infection rate
        Time,...    
        y0);

    y(:,2:4:end) = [];

        
    %%% ODE PART
    
    function y = odemodel(p,nNodes,H,S,F,T,A,W,chi,xi,epsilon,theta,tspan,y0)
        

        [~,y]=ode113(@eqs,tspan,y0);
        
        function dy=eqs(t,y)

            index_t=floor(t-tspan(1))+1;


            dy=zeros(4*nNodes,1);

            P1 = T*y(4:4:end);
           

            P4 = W*(F.*y(4:4:end))./F - sum(W,1)'.*y(4:4:end);

            P4(isnan(P4)) = 0; % for reaches with zero length


            beta_E = p.beta_E;
            theta_C = theta;

            dy(1:4:end) = P1./H.*epsilon - (p.mu_W+p.mu_H)*y(1:4:end);


            dy(2:4:end) = xi.* p.rho_E .* H .* y(1:4:end)./(p.alpha+y(1:4:end))./A ...
                - (p.mu_E + beta_E * S./A).*y(2:4:end);

            
            dy(3:4:end) = beta_E*(1-y(3:4:end)).*y(2:4:end) - p.mu_S*y(3:4:end);


            dy(4:4:end) = theta_C.*y(3:4:end) + P4 - (p.mu_F+chi).*y(4:4:end);

            dy(isnan(dy))=0; % 

   
        end
    end

    
end