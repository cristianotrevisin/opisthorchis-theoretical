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
    s.X = OCN.W;
    out = find(OCN.SC_AccArea == max(OCN.SC_AccArea));
    if OCN.SC_RiverLength(out) == 0
        OCN.SC_RiverLength(out) = OCN.cellsize;
    end
    
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
    % dF is the fish density at the outlet
    ACC_FACTOR = OCN.SC_AccArea.^0.3./max(OCN.SC_AccArea.^0.3);
    s.F = p.dF*OCN.SC_RiverLength.*OCN.SC_RiverWidth.*ACC_FACTOR;
   
    
    % Snail population
    s.S = OCN.SC_RicePaddy_Area*p.dS;
    
   
    % Compute fish mobility matrix
    s.W = calculateWRecursive(1, s.X+s.X', s.F, zeros(size(s.X)), false(s.nNodes,1),1);   
    s.W = s.W*s.par.lambda_F/mean(s.W(:)>0,'all');

    % FISH MARKET
    
    s.chi = s.par.c*s.H;
    s.par.U = sum(s.chi.*s.F,'omitmissing')/sum(s.H,'omitmissing');
    s.delta = -s.chi.*s.F + s.par.U * s.H; s.delta(s.delta<0)=0; % deficit
    s.sigma = s.chi.*s.F - s.par.U * s.H; s.sigma(s.sigma<0) = 0; % surplus
    % 
    % % Build trade matrix
    T = s.delta./exp(OCN.Dist./p.D); T(eye(OCN.nNodes)==1) = 0; % distances
    T = T./sum(T,1).*s.sigma'; % never exceeding surplus
    % 
    % % Check row sum to avoid exceeding deficit import
    T_rows_sum = sum(T,2);
    check_delta = T_rows_sum>s.delta & T_rows_sum>0;
    T(check_delta,:) = T(check_delta,:)./T_rows_sum(check_delta).*s.delta(check_delta);
    clear T_rows_sum
     
    % Add internal fish
    T(eye(OCN.nNodes)==1) = min(s.par.U*s.H, s.chi.*s.F);
     
    s.T = T;
     
    s.A = OCN.SC_RicePaddy_Area; % area for eggs
    s.A(s.A == 0) = NaN;


    % Snail exposure reduction
    s.par.xi = (rand(size(s.H))).*exp(-rand(size(s.H)).*0.00001.*s.H);

    % Raw fish quota
    s.par.epsilon = s.par.epsilon0*ones(s.nNodes,1);

    % Fish exposure
    s.par.theta = s.par.rho_C*s.par.beta_C.*s.S./s.A./(s.par.mu_C+s.par.beta_C.*s.F./s.A);
    s.par.theta(isnan(s.par.theta)) = 0;

    % Unify results
    if ip.Results.Unify == true
        % fish market
        s.chi = s.par.c*sum(s.H.*s.F,'omitmissing')/sum(s.F,'omitmissing'); 
        s.delta = s.chi*sum(s.F,'omitmissing') - s.par.U*sum(s.H,'omitmissing'); s.delta(s.delta<0)=0; % deficit
        s.sigma = -s.chi*sum(s.F,'omitmissing') + s.par.U*sum(s.H,'omitmissing'); s.sigma(s.sigma<0) = 0; % surplus
        s.T = min(s.chi*sum(s.F,'omitmissing'),s.par.U*sum(s.H,'omitmissing')); 

        % populations
        s.nNodes = 1; 
        s.W = 0;
        s.A = sum(s.A,'omitmissing');
        s.H = sum(s.H,'omitmissing');
        s.par.xi = sum(s.par.xi.*s.S,'omitmissing')/sum(s.S,'omitmissing');
        s.par.epsilon = sum(s.par.epsilon.*s.F,'omitmissing')/sum(s.F,'omitmissing');
        s.S = sum(s.S,'omitmissing');
        s.F = sum(s.F,'omitmissing');
        s.par.theta = s.par.rho_C*s.par.beta_C.*s.S/s.A./(s.par.mu_C+s.par.beta_C.*s.F./s.A);
        

        
    end
end