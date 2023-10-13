function [s] = build_setup(OCN,p,TotalPopulation,varargin)

    % Create an inputParser object
    ip = inputParser;

    % Define named properties and their default values
    addRequired(ip, 'OCN');
    addRequired(ip, 'p');
    addRequired(ip, 'TotalPopulation');
    addParameter(ip, 'DownstreamAccumulation', false);
    addParameter(ip, 'Unify', false);
    addParameter(ip, 'seed', nan);

    
 
    % Parse the input arguments
    ip.KeepUnmatched = true;
    parse(ip,OCN,p,TotalPopulation,varargin{:})
    
    if isnan(ip.Results.seed)==false
        rng(ip.Results.seed);
    end

    s.nNodes = OCN.nNodes;

    %%% LOAD PARAMETERS
    s.par = p;
    
    % Hydrological connectivity
    s.W = OCN.W;

    
    %%% POPULATIONS
    % Generate random human population Zipf distributed
    if ip.Results.DownstreamAccumulation == true
        [~,rankacc] = sort(OCN.SC_AccArea,'descend'); ID = zeros(length(rankacc),1);
        for rank = 1:OCN.nNodes
            ID(rankacc(rank)) = rank;
        end
        s.H = zipf(ID', 2, 250);
    else
        s.H = zipf(randperm(OCN.nNodes), 2, 250);
    end
    s.H = s.H/sum(s.H)*TotalPopulation;
    
    % Fish population
    % Scaling 0.3 to flow accumulation, and then multiplying by the catchment's
    % area 
    s.KF = p.dF*OCN.SC_AccArea.^0.3...
        .* OCN.SC_RiverLength.*OCN.SC_RiverWidth.*OCN.SC_RiverDepth...
        /max(OCN.SC_AccArea.^0.3...
        .* OCN.SC_RiverLength.*OCN.SC_RiverWidth.*OCN.SC_RiverDepth)...
        .* OCN.SC_RiverLength.*OCN.SC_RiverWidth.*OCN.SC_RiverDepth;
   
    
    % Snail population
    s.S = OCN.SC_RicePaddy_Area*p.dS;
    
    % FISH MARKET
    %s.par.c = p.U/sum(s.KF); % individual catch rate
    s.par.c = sum(p.U*s.H)/sum(s.H.*s.KF); % individual catch rate
    s.chi = s.par.c*s.H;

    % Fish population at equilibrium
    s.F = find_fish_equilibrium(s.KF,s.W,s.par,s.chi);

    s.delta = 1-s.F*s.par.c/p.U; s.delta(s.delta<0)=0; % deficit
    s.sigma = 1-p.U/s.par.c./s.F; s.sigma(s.sigma<0) = 0; % surplus

    % Build trade matrix
    T = s.delta./exp(OCN.Dist./p.D); T(eye(OCN.nNodes)==1) = 0; % distances
    T = T./sum(T,1).*repmat(s.sigma',OCN.nNodes,1); % never exceeding surplus
    
    % Check row sum to avoid exceeding deficit import
    T_rows_sum = sum(T,2);
    check_delta = T_rows_sum>=s.delta & T_rows_sum>0;
    T(check_delta,:) = T(check_delta,:)./T_rows_sum(check_delta).*s.delta(check_delta);
    clear T_rows_sum
    
    % Add internal fish
    T(eye(OCN.nNodes)==1) = 1-s.sigma;

    s.T = T;

    s.A = OCN.SC_RicePaddy_Area; % area for eggs
    s.A(s.A == 0) = NaN;


    % Add parameter after volume
    betaSH = 9.16e-11;
    s.par.beta_E = betaSH*mean(s.A,'omitmissing')*p.mu_E/(p.rho_E-mean(s.A,'omitmissing')*betaSH*p.dS);

    betaFS = 3.477e-5;
    rhoC = 150;
    MF = 3;
    muC = -log(0.01)/2;
    betaC = betaFS*MF*muC/(rhoC-betaFS*MF*mean(s.F./s.A,'omitmissing'));
    % s.par.theta_C = rhoC*betaC./(muC+betaC*s.F./s.A);
    % s.par.theta_C(isnan(s.par.theta_C)) = 0;
    s.par.theta_C = betaFS*MF;

    % Snail exposure reduction
    s.par.xi = (rand(size(s.H))).*exp(-rand(size(s.H)).*0.00001.*s.H);

    % Unify results
    if ip.Results.Unify == true
        s.chi = s.par.c*sum(s.H.*s.F,'omitmissing')/sum(s.F,'omitmissing'); 
        s.par.c = s.chi/sum(s.H,'omitmissing');
        s.nNodes = 1; 
        s.W = 1;
        s.A = sum(s.A,'omitmissing');
        s.H = sum(s.H,'omitmissing');
        s.KF = sum(s.KF,'omitmissing');
        s.par.xi = sum(s.par.xi.*s.S,'omitmissing')/sum(s.S,'omitmissing');
        s.S = sum(s.S,'omitmissing');
        s.F = sum(s.F,'omitmissing');
        s.par.theta_C = rhoC*betaC./(muC+betaC*s.F./s.A);
        s.delta = 1-s.F*s.par.c/p.U; s.delta(s.delta<0)=0; % deficit
        s.sigma = 1-p.U/s.par.c./s.F; s.sigma(s.sigma<0) = 0; % surplus
        s.T = diag(1-s.sigma); 
    end
end