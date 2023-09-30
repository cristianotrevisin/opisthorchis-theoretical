function par = common_parameters()

    % Define parameters
    par.alpha_M = 1;
    par.mu_H = 1/(365.25*68.5); % Life expectancy Laos
    par.mu_W = 1/(365.25*10); % Burli
    par.gamma = 1/(365.25*2); % FIND
    par.alpha_E = 0.1; % amount of eggs that are dumped in snail habitat %%% NEW!!!!!!
    par.rho_E = par.alpha_E*1e6; % Crellen (10000 for a 5000 worms infection)
    par.omega = 100; % Crellen (10 for 1 worm infection)
    par.mu_E = -log(0.8)/21; % 
    par.mu_S = 1/(365.25); %  BURLI 
    par.theta_C = 1e-9; % FIND
    par.mu_F = 1/(365.25*2.5); % BURLI
    par.U = 1;
    par.D = 10000;
    
    par.dF = 10;
    par.dS = 30;
    par.lambda_FU=0;
    par.lambda_FD=0;


end