function [UP,DOWN,fields] = sensitivity_analysis_EE(OCN,p,seed)
    
    fields = fieldnames(p);
    UP = zeros(4,length(fields));
    DOWN = zeros(4,length(fields));

    % Run baseline
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    bsl = get_equilibrium_values(s);


    for pp = 1:length(fields)
        pp
        p.(fields{pp}) = p.(fields{pp})*1.2;
        s = build_setup(OCN,p,33*800*1000,'seed',seed);
        tmp = get_equilibrium_values(s);
        UP(:,pp) = (tmp-bsl)./bsl;
        
        p.(fields{pp}) = p.(fields{pp})/1.2*0.8;
        s = build_setup(OCN,p,33*800*1000,'seed',seed);
        tmp = get_equilibrium_values(s);
        DOWN(:,pp) = (tmp-bsl)./bsl;
        p.(fields{pp}) = p.(fields{pp})/0.8;
    end

   

    % xi
    s.par.xi = s.par.xi * 1.2;
    s.par.xi(s.par.xi>1) = 1;
    tmp = get_equilibrium_values(s);
    UP(:,end+1) = (tmp-bsl)./bsl;

    s.par.xi = s.par.xi / 1.2*0.8;
    tmp = get_equilibrium_values(s);
    DOWN(:,end+1) = (tmp-bsl)./bsl;
    s.par.xi = s.par.xi / 0.8;
    fields{end+1} = '\xi';


    taboo = ["D", "lambda_F"];
    for i = 1:length(taboo)
        tmp = find(strcmp(fields,taboo(i)));
        if ~isempty(tmp)
            UP(:,tmp) = [];
            DOWN(:,tmp) = [];
            fields(tmp) = [];
        end
    end



    % Checking what happens with distance D
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    seeding_node = find(s.sigma == max(s.sigma));
    bsl = get_simulation_equilibrium(s,seeding_node);

    p.D = p.D*1.2;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    tmp = get_simulation_equilibrium(s,seeding_node);
    UP(:,end+1) = (tmp-bsl)./bsl;

    p.D = p.D/1.2*0.8;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    tmp = get_simulation_equilibrium(s,seeding_node);
    DOWN(:,end+1) = (tmp-bsl)./bsl;
    fields{end+1} = 'D';
    p.D = p.D/0.8;

    % Checking what happens with parameters lambda
    outlet = find(OCN.SC_AccArea == max(OCN.SC_AccArea));
    seeding_node = find(OCN.distW(outlet,:)==max(OCN.distW(outlet,:)));
    bsl = get_simulation_equilibrium(s,seeding_node);

    p.lambda_F = p.lambda_F*1.2;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    tmp = get_simulation_equilibrium(s,seeding_node);
    UP(:,end+1) = (tmp-bsl)./bsl;

    p.lambda_F = p.lambda_F/1.2*0.8;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    tmp = get_simulation_equilibrium(s,seeding_node);
    DOWN(:,end+1) = (tmp-bsl)./bsl;
    p.lambda_F = p.lambda_F/0.8;
    fields{end+1} = 'lambda_{F}';

 

    % Format labels for easier plotting
    fields = strrep(fields,'mu','\mu');
    fields = strrep(fields,'gamma','\gamma');
    fields = strrep(fields,'rho','\rho');
    fields = strrep(fields,'omega','\omega');
    fields = strrep(fields,'theta','\theta');
    fields = strrep(fields,'lambda','\lambda');
    fields = strrep(fields,'beta','\beta');
    fields = strrep(fields,'dF','d_F');
    fields = strrep(fields,'dS','d_S');
    fields = strrep(fields,'alpha','\alpha');
    fields = strrep(fields,'epsilon0','\epsilon_0');
    fields = strrep(fields,'c','k');

end

function out = get_equilibrium_values(setup)
    setup.A(isnan(setup.A)) = 0; setup.S(isnan(setup.S)) = 0; setup.F(isnan(setup.F)) = 0;
   
    eqi = zeros(setup.nNodes,4);
    for nn = 1:setup.nNodes
        tem = max(find_EE(setup.par,setup.H(nn), setup.S(nn), setup.A(nn),...
            setup.chi(nn),setup.par.epsilon(nn),setup.par.theta(nn),...
            setup.par.xi(nn),setup.T(nn,nn)));
        if  isempty(tem) || tem(1)==0
            eqi(nn,:)=0;
        else
            eqi(nn,:)=tem;
        end
    end

    out = [eqi(:,1)'*setup.H/sum(setup.H); 
        eqi(:,2)'*setup.A/sum(setup.A); 
        eqi(:,3)'*setup.S/sum(setup.S);
        eqi(:,4)'*setup.F/sum(setup.F)];

end

function out = get_simulation_equilibrium(setup,seeding_node)
    Time = 1:1000*365;
    setup.period = ones(length(Time),1);
    
    y0 = zeros(setup.nNodes,4);

    y0(seeding_node,1) = 100/setup.H(seeding_node);

    
    y = model_ODE(Time,setup.par,setup,y0');
    setup.A(isnan(setup.A)) = 0; setup.S(isnan(setup.S)) = 0; setup.F(isnan(setup.F)) = 0;
    out = [(y(end,1:4:end)*setup.H)/sum(setup.H); 
        (y(end,2:4:end)*setup.A)/sum(setup.A);...
        (y(end,3:4:end)*setup.S)/sum(setup.S);
        (y(end,4:4:end)*setup.F)/sum(setup.F)];

end