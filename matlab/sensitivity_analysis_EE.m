function [UP,DOWN,fields] = sensitivity_analysis_EE(OCN,p,seed)
    
    fields = fieldnames(p);
    UP = zeros(3,length(fields),3,2);
    DOWN = zeros(3,length(fields));

    f_up = [1.2 1.1 1.05];
    f_down = [0.8 0.9 0.95];

    % Run baseline
    p.lambda_F = 0;
    s1 = build_setup(OCN,p,33*800*1000,'seed',seed);
    s1.T = diag(diag(s1.T));
    bsl1 = get_simulation_equilibrium(s1);

    s2 = build_setup(OCN,p,33*800*1000,'seed',seed,'DownstreamAccumulation',true);
    s2.T = diag(diag(s2.T));
    bsl2 = get_simulation_equilibrium(s2);


    for pp = 1:length(fields)-2
        pp
        
        for fm = 1:3
            p.(fields{pp}) = p.(fields{pp})*f_up(fm);
            s1 = build_setup(OCN,p,33*800*1000,'seed',seed);
            s1.T = diag(diag(s1.T));
            tmp = get_simulation_equilibrium(s1);
            UP(:,pp,fm,1) = (tmp-bsl1)./bsl1;
            s2 = build_setup(OCN,p,33*800*1000,'seed',seed,'DownstreamAccumulation',true);
            s2.T = diag(diag(s2.T));
            tmp = get_simulation_equilibrium(s2);
            UP(:,pp,fm,2) = (tmp-bsl2)./bsl2;
            
            p.(fields{pp}) = p.(fields{pp})/f_up(fm)*f_down(fm);
            s1 = build_setup(OCN,p,33*800*1000,'seed',seed);
            s1.T = diag(diag(s1.T));
            tmp = get_simulation_equilibrium(s1);
            DOWN(:,pp,fm,1) = (tmp-bsl1)./bsl1;
            s2 = build_setup(OCN,p,33*800*1000,'seed',seed,'DownstreamAccumulation',true);
            s2.T = diag(diag(s2.T));
            tmp = get_simulation_equilibrium(s2);
            DOWN(:,pp,fm,2) = (tmp-bsl2)./bsl2;
            p.(fields{pp}) = p.(fields{pp})/f_down(fm);
        end
    end

    s1 = build_setup(OCN,p,33*800*1000,'seed',seed);
    s1.T = diag(diag(s1.T));
    s2 = build_setup(OCN,p,33*800*1000,'seed',seed,'DownstreamAccumulation',true);
    s2.T = diag(diag(s2.T));
    xi_1 = s1.par.xi;
    xi_2 = s2.par.xi;

    LF = length(fields);

    % xi
    for fm = 1:3
        s1.par.xi = xi_1 * f_up(fm);
        s1.par.xi(s1.par.xi>1) = 1;
        tmp = get_simulation_equilibrium(s1);
        UP(:,LF+1,fm,1) = (tmp-bsl1)./bsl1;
        s2.par.xi = xi_2 * f_up(fm);
        s2.par.xi(s2.par.xi>1) = 1;
        tmp = get_simulation_equilibrium(s2);
        UP(:,LF+1,fm,2) = (tmp-bsl2)./bsl2;
    
        s1.par.xi = xi_1 * f_down(fm);
        tmp = get_simulation_equilibrium(s1);
        DOWN(:,LF+1,fm,1) = (tmp-bsl1)./bsl1;
        s2.par.xi = xi_2 * f_down(fm);
        tmp = get_simulation_equilibrium(s2);
        DOWN(:,LF+1,fm,2) = (tmp-bsl2)./bsl2;
        fields{LF+1} = '\xi';
    end

 
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


function out = get_simulation_equilibrium(setup)
    Time = 1:1000*365;
    
    y0 = zeros(setup.nNodes,4);

    y0(:,1) = 0.1;

    
    y = model_ODE(Time,setup.par,setup,y0');
    setup.S(isnan(setup.S)) = 0; setup.F(isnan(setup.F)) = 0;
    out = [(y(end,1:3:end)*setup.H)/sum(setup.H); 
        (y(end,2:3:end)*setup.S)/sum(setup.S);
        (y(end,3:3:end)*setup.F)/sum(setup.F)];

end