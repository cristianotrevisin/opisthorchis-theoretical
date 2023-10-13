function y = model_ODE(Time,par,setup,y0)

    %node_out = find(sum(W,1)==0);
    outlet = zeros(setup.nNodes,1); %outlet(sum(setup.W,1)==0)=1;
    y=odemodel(par,...
        setup.nNodes,... #number of nodes
        setup.H,... #human population
        setup.S,... #snail population
        setup.F,...#carrying fish population
        setup.T,... #caught fish trade
        setup.A,... #local area for snails
        setup.W,... #hydrological connectivity
        setup.chi,...#fish catch rate
        outlet,...#allowing release outside of the system
        setup.par.xi,... #urbanization reduction of snails' exposure
        setup.period,...#periodic function
        Time,...    
        y0);
        
    %%% ODE PART
    
    function y = odemodel(p,nNodes,H,S,F,T,A,W,chi,outlet,xi,period,tspan,y0)

        [~,y]=ode45(@eqs,tspan,y0);
        
        function dy=eqs(t,y)

            index_t=floor(t-tspan(1))+1;


            dy=zeros(4*nNodes,1);

            P1 = T*(F.*y(4:4:end).*chi);
            

            P4 = p.lambda_FU*W'*(F.*y(4:4:end))./F ...
               + p.lambda_FD*W*(F.*y(4:4:end))./F ...
               - (p.lambda_FD*(sum(W,1)'+outlet) + p.lambda_FU*sum(W,2)).*y(4:4:end);

            P4(isnan(P4)) = 0; % for reaches with zero length


            beta_E = p.beta_E*period(index_t);
            theta_C = p.theta_C*period(index_t);

            dy(1:4:end) = P1./H - (p.mu_H+p.mu_W+p.gamma)*y(1:4:end);


            dy(2:4:end) = xi.* p.rho_E .* H .* y(1:4:end)./(p.omega+y(1:4:end))./A ...
                - (p.mu_E + beta_E * S./A).*y(2:4:end);

            
            dy(3:4:end) = beta_E*(1-y(3:4:end)).*y(2:4:end) - p.mu_S*y(3:4:end);


            dy(4:4:end) = theta_C.*S.*y(3:4:end) + P4 - (p.mu_F+chi).*y(4:4:end);

            dy(isnan(dy))=0; % because the volume is 0 in a few instances

   
        end
    end

    
end