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
    %p.c = sum(p.U*s.H)/sum(s.H.*s.KF); % individual catch rate
    s.delta = 1-s.KF*p.c/p.U; s.delta(s.delta<0)=0; % deficit
    s.sigma = 1-p.U/p.c./s.KF; s.sigma(s.sigma<0) = 0; % surplus
    s.chi = p.c*s.H;
    
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

    s.V = OCN.SC_RicePaddy_Area.*0.5; % volume for eggs

    % Hydrological connectivity
    s.W = OCN.W;
    


    % Unify results
    if ip.Results.Unify == true
        s.chi = p.c*sum(s.H); 
        s.nNodes = 1; 
        s.T = 1; 
        s.W = 1;
        s.V = sum(s.V);
        s.H = sum(s.H);
        s.KF = sum(s.KF);
        s.S = sum(s.S);
    end
end