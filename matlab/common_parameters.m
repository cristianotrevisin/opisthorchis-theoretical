function par = common_parameters()

    % Define parameters
    % Mortality
    par.mu_H = 1/(365.25*68.5); % Life expectancy Laos
    par.mu_W = 1/(365.25*10); % Burli
    par.mu_E = -log(0.8)/21; % 
    par.mu_S = 1/(365.25); %  BURLI 
    par.mu_F = 1/(365.25*2.5); % BURLI

    %Recovery
    par.gamma = 1/(365.25*2); % FIND

    %Egg shedding
    par.rho_E = 102.4353*200; % Crellen 
    par.omega = 424.7445; % Crellen (10 for 1 worm infection)

    par.U = .1;
    par.D = 10000;
    
    par.dF = 0.08;
    par.dS = 20;
    par.lambda_FU=0.01;
    par.lambda_FD=0.01;


end