function [I,S,F] = EE_OPI(par,setup,n)

    betaH = par.beta_HS*setup.H(n);
    betaS = par.beta_SF*setup.Ns(n);
    betaF = par.beta_FH*setup.Nf(n)*setup.Cf(n);

    muH = par.mu_H;
    muS = par.mu_S;
    muF = par.mu_F + setup.Cf(n);

    if betaH*betaS*betaF-muH*muS*muF <0
        raise('Attention, Endemic Equilibrium does not exist')
    end

    I = (betaH*betaS*betaF-muH*muS*muF)/betaH/(betaF*betaS+betaS*muH+muH*muF);
    S = (betaH*betaS*betaF-muH*muS*muF)/betaS/(betaF*betaH+betaF*muS+muH*muS);
    F = (betaH*betaS*betaF-muH*muS*muF)/betaF/(betaH*betaS+betaH*muF+muS*muF);

end