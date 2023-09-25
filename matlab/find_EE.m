function sol = find_EE(p,chi,H,S,KF,V)

    IH = sym('IH','real');
    E = sym('E','real');
    YS = sym('YS','real');
    F = sym('F','real');
    IF = sym('IF','real');

    eq1 = p.alpha_M*chi*F/H*IF - (p.mu_H+p.mu_W+p.gamma)*IH==0;   



    eq2 = p.rho_E/V*H*IH/(p.omega+IH) - (p.mu_E+p.beta_E*S/V)*E == 0;  

    eq3 = p.beta_E*E.*(1-YS) - p.mu_S.*YS == 0; 

    eq4 = p.mu_F*(KF-F) - chi*F == 0;

    eq5 = p.theta_C*S*YS - (p.mu_F+chi)*IF;

    sol = solve([eq1;eq2;eq3;eq4;eq5],[IH,E,YS,F,IF]);

    sol = [double(sol.IH) double(sol.E) double(sol.YS) ...
        double(sol.F) double(sol.IF)];


    % sol = zeros(4,1);
    % NUM = p.alpha_M*chi*p.rho_E/V*p.beta_E*p.theta_C*S - ...
    %     p.omega*(p.mu_H+p.gamma+p.mu_W)*(p.mu_E+p.beta_E*S/V)*...
    %     p.mu_S*(p.mu_F+chi)
    % 
    % DEN1 = (p.mu_H+p.gamma+p.mu_W)*(p.mu_F+chi)*...
    %     (p.rho_E/V*H*p.beta_E+p.mu_S*(p.mu_E+p.beta_E*S/V));
    % 
    % DEN2 = p.beta_E*(p.mu_E+p.beta_E*S/V)*(p.alpha_M*chi*F/H*p.theta_C*S +...
    %     p.omega*(p.mu_H+p.gamma+p.mu_W)*(p.mu_F+chi));
    % 
    % DEN3 = p.alpha_M*chi*F/H*p.theta_C*S*(p.rho_E/V*H*p.beta_E+...
    %     (p.mu_E+p.beta_E*S/V)*p.mu_S);
    % 
    % DEN4 = p.alpha_M*chi*F/H*(p.mu_F+chi)*(p.rho_E/V*p.beta_E+p.mu_S*...
    %     (p.mu_E+p.beta_E*S/V));
    % 
    % if NUM > 0
    %     sol(1) = NUM/DEN1;
    %     sol(2) = NUM/DEN2;
    %     sol(3) = NUM/DEN3;
    %     sol(4) = NUM/DEN4;
    % else
    %     sol(1:end) = NaN;
    % end


end