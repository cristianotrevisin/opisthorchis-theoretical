function F = find_fish_equilibrium(KF,W,par,chi)

    FS = sym('FS',size(KF),'real');

    eqs = par.mu_F*(KF-FS)-chi.*FS + ...
        par.lambda_FU*W'*FS + par.lambda_FD*W*FS ...
        - (par.lambda_FD*(sum(W,1)') + par.lambda_FU*sum(W,2)).*FS == 0;

    sol = struct2cell(solve(eqs,FS));

    F = double([sol{:}])';

end