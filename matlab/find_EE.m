function sol = find_EE(par,setup)

    WH = sym('WH',[setup.nNodes 1]);
    IS = sym('IS',[setup.nNodes 1]);
    IF = sym('IF',[setup.nNodes 1]);
    
    eq1 = par.beta_FH*(setup.M*(IF.*setup.Nf)) - par.mu_H*WH==0;     

    eq2 = par.beta_HS*setup.H.*WH.*(1-IF) - par.mu_S*IS == 0;  

    eq3 = par.beta_SF.*setup.Ns.*IS.*(1-IF) - par.mu_F.*IF == 0; 

    sol = solve([eq1;eq2;eq3],[WH,IS,IF]);

end