function sol = compute_EE(s)
% FIND_EE  Find the equilibria.
    p = s.par;

    a1 = s.T.*p.epsilon./s.H;               a1(isnan(a1))=0;
    a2 = p.xi*p.rho_E./s.A.*s.H;            a2(isnan(a2))=0;
    a3 = p.beta_E;                          a3(isnan(a3))=0;
    a4 = p.theta;                           a4(isnan(a4))=0;

    b1 = p.mu_H+p.mu_W;             b1(isnan(b1))=0;
    b2 = p.mu_E+p.beta_E.*s.S./s.A;         b2(isnan(b2))=0;
    b3 = p.mu_S;                            b3(isnan(b3))=0;
    b4 = p.mu_F+s.chi;                      b4(isnan(b4))=0;

    NUM = a1.*a2.*a3.*a4 - p.alpha.*b1.*b2.*b3.*b4;

    NUM(NUM<0)=0;


    DEN = a2.*a3.*b1.*b4+b1.*b2.*b3.*b4;
   
    X1 = NUM./DEN;
    X4 = X1.*b1./a1;
    X3 = X4.*b4./a4;

    sol = [X1'; X3'; X4']';

    sol(isnan(sol))=0;





end