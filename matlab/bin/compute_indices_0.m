function [R, E, T, S] = compute_indices(par,setup,y,Time)

% unwrap
WH = y(:,1:3:end);
SI = y(:,2:3:end);
FI = y(:,3:3:end);

R = zeros(1,length(Time));
E = zeros(1,length(Time));


% Necessary vectors
n=setup.nNodes;
u=ones(1,n); 
u1n=1:n; 
U=sparse(u1n,u1n,u,n,n); 
Z=sparse(zeros(n)); 

for t = 1:200:Time(end)
    t
    T = zeros(setup.nNodes*3,setup.nNodes*3);
    S = zeros(setup.nNodes*3,setup.nNodes*3);
    %Build transmission matrix
    for n = 1:setup.nNodes
        % human equations
        T(1+3*(n-1),1+3*(n-1)) = -par.beta_FH*sum(setup.M(n,:).*setup.Nf'.*setup.Cf'.*FI(t,:)); 
        T(1+3*(n-1),3+3*(n-1)) =  par.beta_FH*setup.Nf(n)*setup.Cf(n)*(1-WH(t,n)); 
        S(1+3*(n-1),1+3*(n-1)) = -par.mu_H;

        % snails equations
        S(2+3*(n-1),2+3*(n-1)) = -par.beta_HS*setup.H(n)*WH(t,n); 
        S(2+3*(n-1),1+3*(n-1)) =  par.beta_HS*setup.H(n)*(1-SI(t,n)); 
        S(2+3*(n-1),2+3*(n-1)) = -par.mu_S - par.mS;

        % fish equations
        S(3+3*(n-1),3+3*(n-1)) = -par.beta_SF*setup.Ns(n)*SI(t,n); 
        S(3+3*(n-1),2+3*(n-1)) =  par.beta_SF*setup.Ns(n)*(1-FI(t,n)); 
        S(3+3*(n-1),3+3*(n-1)) = -par.mu_F - par.mF - setup.Cf(n);

        for m = 1:setup.nNodes
%         human equations
        T(1+3*(n-1),3+3*(m-1)) = par.beta_FH*setup.M(n,m)*setup.Nf(m)*setup.Cf(m)*(1-WH(t,n));

%         snails equations
        S(2+3*(n-1),2+3*(m-1)) = S(2+3*(n-1),2+3*(m-1))+...
            par.mS*setup.W(n,m)*setup.Ns(m); 

%         fish equations
        S(3+3*(n-1),3+3*(m-1)) = S(3+3*(n-1),3+3*(m-1))+...
            par.mF*setup.W2(n,m)*setup.Nf(m); 
        end
    end

    J = T+S;
    H = (J+J')/2;
    NGM = -T*inv(S);


    R(t) = max(real(eig(full(NGM))));
    E(t) = eigs(H,1,'largestreal');

end

figure()
subplot(3,1,1)
plot(Time,R)
ylabel('$\mathcal{R}_t$','Interpreter','latex')
subplot(3,1,2)
plot(Time,E)
ylabel('$e_t$','Interpreter','latex')
subplot(3,1,3)
hold on
plot(Time,WH*setup.H/sum(setup.H))
plot(Time,SI*setup.Ns/sum(setup.Ns))
plot(Time,FI*setup.Nf/sum(setup.Nf))
legend('I','S','F')
xlabel('Time')


end