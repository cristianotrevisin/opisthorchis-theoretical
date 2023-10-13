function [UP,DOWN,fields] = sensitivity_analysis_EE(OCN,p,seed)
    
    fields = fieldnames(p);
    UP = zeros(4,length(fields));
    DOWN = zeros(4,length(fields));

    % Run baseline
    s = build_setup(OCN,p,33*800*1000,'seed',seed,'Unify',true);
    bsl = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';

    for pp = 1:length(fields)
        p.(fields{pp}) = p.(fields{pp})*1.2;
        s = build_setup(OCN,p,33*800*1000,'seed',seed,'Unify',true);
        tmp = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';
        UP(:,pp) = (tmp-bsl)./bsl;
        
        p.(fields{pp}) = p.(fields{pp})/1.2*0.8;
        s = build_setup(OCN,p,33*800*1000,'seed',seed,'Unify',true);
        tmp = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';
        DOWN(:,pp) = (tmp-bsl)./bsl;
        p.(fields{pp}) = p.(fields{pp})/0.8;
    end

    %Other parameters that should be checked: beta_E, theta_C, c, xi
    % beta_E
    s.par.beta_E = s.par.beta_E * 1.2;
    tmp = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';
    UP(:,end+1) = (tmp-bsl)./bsl;
    s.par.beta_E = s.par.beta_E / 1.2*0.8;
    tmp = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';
    DOWN(:,end+1) = (tmp-bsl)./bsl;
    s.par.beta_E = s.par.beta_E / 0.8;
    fields{end+1} = 'beta_E';

    % theta_C
    s.par.theta_C = s.par.theta_C * 1.2;
    tmp = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';
    UP(:,end+1) = (tmp-bsl)./bsl;
    s.par.theta_C = s.par.theta_C / 1.2*0.8;
    tmp = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';
    DOWN(:,end+1) = (tmp-bsl)./bsl;
    s.par.theta_C = s.par.theta_C / 0.8;
    fields{end+1} = 'theta_C';

    % c
    s.chi = s.chi*1.2;
    tmp = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';
    UP(:,end+1) = (tmp-bsl)./bsl;
    s.chi = s.chi / 1.2*0.8;
    tmp = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';
    DOWN(:,end+1) = (tmp-bsl)./bsl;
    s.chi = s.chi / 0.8;
    fields{end+1} = '\chi';

    % xi
    s.par.xi = s.par.xi * 1.2;
    tmp = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';
    UP(:,end+1) = (tmp-bsl)./bsl;
    s.par.xi = s.par.xi / 1.2*0.8;
    tmp = max(find_EE(s.par,s.chi,s.H,s.S,s.F,s.A,s.par.theta_C,s.par.xi))';
    DOWN(:,end+1) = (tmp-bsl)./bsl;
    s.par.xi = s.par.xi / 0.8;
    fields{end+1} = '\xi';


    taboo = ["D", "lambda_FD", "lambda_FU"];
    for i = 1:length(taboo)
        tmp = find(strcmp(fields,taboo(i)));
        if ~isempty(tmp)
            UP(:,tmp) = [];
            DOWN(:,tmp) = [];
            fields(tmp) = [];
        end
    end



    % Checking what happens with distance D
    Time = 1:40*365;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    s.period = ones(length(Time),1);
    surplus = s.sigma.*s.F;
    
    y0 = zeros(OCN.nNodes,4);
    y0(surplus==max(surplus),1) = 1/s.H(surplus==max(surplus));
    
    y = model_ODE(Time,s.par,s,y0');
         s.A(isnan(s.A)) = 0; s.S(isnan(s.S)) = 0; s.F(isnan(s.F)) = 0;
    bsl = [max(y(:,1:4:end)*s.H)/sum(s.H); max(y(:,2:4:end)*s.A)/sum(s.A);...
        max(y(:,3:4:end)*s.S)/sum(s.S); max(y(:,4:4:end)*s.F)/sum(s.F)];

    p.D = p.D*1.2;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    s.period = ones(length(Time),1);
    y = model_ODE(Time,s.par,s,y0');
         s.A(isnan(s.A)) = 0; s.S(isnan(s.S)) = 0; s.F(isnan(s.F)) = 0;
    tmp = [max(y(:,1:4:end)*s.H)/sum(s.H); max(y(:,2:4:end)*s.A)/sum(s.A);...
        max(y(:,3:4:end)*s.S)/sum(s.S); max(y(:,4:4:end)*s.F)/sum(s.F)];
    UP(:,end+1) = (tmp-bsl)./bsl;

    p.D = p.D/1.2*0.8;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    s.period = ones(length(Time),1);
    y = model_ODE(Time,s.par,s,y0');
         s.A(isnan(s.A)) = 0; s.S(isnan(s.S)) = 0; s.F(isnan(s.F)) = 0;
    tmp = [max(y(:,1:4:end)*s.H)/sum(s.H); max(y(:,2:4:end)*s.A)/sum(s.A);...
        max(y(:,3:4:end)*s.S)/sum(s.S); max(y(:,4:4:end)*s.F)/sum(s.F)];
    DOWN(:,end+1) = (tmp-bsl)./bsl;
    fields{end+1} = 'D';
    p.D = p.D/0.8;

    % Checking what happens with parameters lambda
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    s.period = ones(length(Time),1);
    outlet = find(OCN.SC_AccArea == max(OCN.SC_AccArea));
    seeding_node = find(OCN.distW(outlet,:)==max(OCN.distW(outlet,:)));
    y0 = zeros(OCN.nNodes,4);
    y0(seeding_node,1) = 1/s.H(seeding_node);
    y = model_ODE(Time,s.par,s,y0');
         s.A(isnan(s.A)) = 0; s.S(isnan(s.S)) = 0; s.F(isnan(s.F)) = 0;
    bsl = [max(y(:,1:4:end)*s.H)/sum(s.H); max(y(:,2:4:end)*s.A)/sum(s.A);...
        max(y(:,3:4:end)*s.S)/sum(s.S); max(y(:,4:4:end)*s.F)/sum(s.F)];

    p.lambda_FD = p.lambda_FD*1.2;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    s.period = ones(length(Time),1);
    y = model_ODE(Time,s.par,s,y0');
         s.A(isnan(s.A)) = 0; s.S(isnan(s.S)) = 0; s.F(isnan(s.F)) = 0;
    tmp = [max(y(:,1:4:end)*s.H)/sum(s.H); max(y(:,2:4:end)*s.A)/sum(s.A);...
        max(y(:,3:4:end)*s.S)/sum(s.S); max(y(:,4:4:end)*s.F)/sum(s.F)];
    UP(:,end+1) = (tmp-bsl)./bsl;

    p.lambda_FD = p.lambda_FD/1.2*0.8;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    s.period = ones(length(Time),1);
    y = model_ODE(Time,s.par,s,y0');
         s.A(isnan(s.A)) = 0; s.S(isnan(s.S)) = 0; s.F(isnan(s.F)) = 0;
    tmp = [max(y(:,1:4:end)*s.H)/sum(s.H); max(y(:,2:4:end)*s.A)/sum(s.A);...
        max(y(:,3:4:end)*s.S)/sum(s.S); max(y(:,4:4:end)*s.F)/sum(s.F)];
    DOWN(:,end+1) = (tmp-bsl)./bsl;
    p.lambda_FD = p.lambda_FD/0.8;

    fields{end+1} = 'lambda_{FD}';

    p.lambda_FU = p.lambda_FU*1.2;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    s.period = ones(length(Time),1);
    y = model_ODE(Time,s.par,s,y0');
         s.A(isnan(s.A)) = 0; s.S(isnan(s.S)) = 0; s.F(isnan(s.F)) = 0;
    tmp = [max(y(:,1:4:end)*s.H)/sum(s.H); max(y(:,2:4:end)*s.A)/sum(s.A);...
        max(y(:,3:4:end)*s.S)/sum(s.S); max(y(:,4:4:end)*s.F)/sum(s.F)];
    UP(:,end+1) = (tmp-bsl)./bsl;

    p.lambda_FU = p.lambda_FU/1.2*0.8;
    s = build_setup(OCN,p,33*800*1000,'seed',seed);
    s.period = ones(length(Time),1);
    y = model_ODE(Time,s.par,s,y0');
         s.A(isnan(s.A)) = 0; s.S(isnan(s.S)) = 0; s.F(isnan(s.F)) = 0;
    tmp = [max(y(:,1:4:end)*s.H)/sum(s.H); max(y(:,2:4:end)*s.A)/sum(s.A);...
        max(y(:,3:4:end)*s.S)/sum(s.S); max(y(:,4:4:end)*s.F)/sum(s.F)];
    DOWN(:,end+1) = (tmp-bsl)./bsl;
    p.lambda_FU = p.lambda_FU/0.8;

    fields{end+1} = 'lambda_{FU}';

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





   
end