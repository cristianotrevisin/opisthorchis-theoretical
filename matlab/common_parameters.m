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
    par.alpha = 424.7445; % Crellen (10 for 1 worm infection)

    par.D = 10000;
    
    par.rho_C = (263*6+98*9+91*17)/(6+9+17); % cercariae/day
    par.mu_C = -log(0.01)/2; 


    par.dF = .08;
    par.epsilon0 = 0.5;

    par.beta_E =  .7e-07;
    par.beta_C = .41;
    par.c = 1e-7;


    par.dS = 20;
    par.lambda_F=0;


end