function sol = find_EE(p,chi,H,S,F,V)
% FIND_EE  Find the equilibria.

    IH = sym('IH','real');
    E = sym('E','real');
    YS = sym('YS','real');
    IF = sym('IF','real');

    eq1 = p.alpha_M*chi*F/H*IF - (p.mu_H+p.mu_W+p.gamma)*IH==0;   

    eq2 = p.rho_E/V*H*IH/(p.omega+IH) - (p.mu_E+p.beta_E*S/V)*E == 0;  

    eq3 = p.beta_E*E.*(1-YS) - p.mu_S.*YS == 0; 

    eq4 = p.theta_C*S*YS - (p.mu_F+chi)*IF;

    sol = solve([eq1;eq2;eq3;eq4],[IH,E,YS,IF]);

    sol = [double(sol.IH) double(sol.E) double(sol.YS) double(sol.IF)];

end