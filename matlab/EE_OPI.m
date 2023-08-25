function [I,S,F] = EE_OPI(par,setup,n)

    betaH = par.beta_HS*setup.H(n);
    betaS = par.beta_SF*setup.Ns(n);
    betaF = par.beta_FH*setup.Nf(n);

    muH = par.mu_H;
    muS = par.mu_S;
    muF = par.mu_F;
betaH*betaS*betaF-muH*muS*muF
    if betaH*betaS*betaF-muH*muS*muF <0
        warning('Attention, Endemic Equilibrium does not exist')
        eene = true;
    else
        eene = false;
    end

    I = (betaH*betaS*betaF-muH*muS*muF)/betaH/muH/(betaS+muF);
    S = (betaH*betaS*betaF-muH*muS*muF)/betaS/(betaF*betaH+muH*muS);
    F = (betaH*betaS*betaF-muH*muS*muF)/betaF/betaH/(betaS+muF);
    
    if eene
        I = NaN; S = NaN; F = NaN;
    end

end