function y = model_ODE(Time,par,setup,y0)

    %node_out = find(sum(W,1)==0);

    y=odemodel(par,...
        setup.nNodes,... #number of nodes
        setup.H,... #human population
        setup.S,... #snail population
        setup.F,...#carrying fish population
        setup.T,... #caught fish trade
        setup.A,... #local area for snails
        setup.W,... #hydrological connectivity
        setup.chi,...#fish catch rate
        setup.par.xi,... #urbanization reduction of snails' exposure
        setup.par.epsilon,... #fraction of fish consumed raw
        setup.par.theta,... #fish infection rate
        Time,...    
        y0);

    y(:,2:4:end) = [];

        
    %%% ODE PART
    
    function y = odemodel(p,nNodes,H,S,F,T,A,W,chi,xi,epsilon,theta,tspan,y0)
        

        [~,y]=ode45(@eqs,tspan,y0);
        
        function dy=eqs(t,y)

            index_t=floor(t-tspan(1))+1;


            dy=zeros(4*nNodes,1);

            P1 = T*y(4:4:end);
           

            P4 = W*(F.*y(4:4:end))./F - sum(W,1)'.*y(4:4:end);

            P4(isnan(P4)) = 0; % for reaches with zero length


            beta_E = p.beta_E;
            theta_C = theta;

            dy(1:4:end) = P1./H.*epsilon - (p.mu_W+p.mu_H)*y(1:4:end);


            dy(2:4:end) = xi.* p.rho_E .* H .* y(1:4:end)./(p.alpha+y(1:4:end))./A ...
                - (p.mu_E + beta_E * S./A).*y(2:4:end);

            
            dy(3:4:end) = beta_E*(1-y(3:4:end)).*y(2:4:end) - p.mu_S*y(3:4:end);


            dy(4:4:end) = theta_C.*y(3:4:end) + P4 - (p.mu_F+chi).*y(4:4:end);

            dy(isnan(dy))=0; % 

   
        end
    end

    
end