function y = model_ODE(Time,par,setup,y0)

    %node_out = find(sum(W,1)==0);
    outlet = zeros(setup.nNodes,1); %outlet(sum(setup.W,1)==0)=1;
    y=odemodel(par,...
        setup.nNodes,... #number of nodes
        setup.H,... #human population
        setup.S,... #snail population
        setup.KF,...#carrying fish population
        setup.T,... #caught fish trade
        setup.V,... #local volumes
        setup.W,... #hydrological connectivity
        setup.chi,...#fish catch rate
        outlet,...#allowing release outside of the system
        Time,...    
        y0);
        
    %%% ODE PART
    
    function y = odemodel(p,nNodes,H,S,KF,T,V,W,chi,outlet,tspan,y0)

        [~,y]=ode45(@eqs,tspan,y0);
        
        function dy=eqs(t,y)

            % NOTE -> We assume a constant trade matrix to reduce
            % computational effort
            
            dy=zeros(5*nNodes,1);

            P1 = T*(y(4:5:end).*y(5:5:end).*chi);
            

            if any(isnan(P1))
                warning('There are nans')
            end

            %P2 = (W*(y(2:5:end).*V))./V - sum(W,1)'.*y(2:5:end);
        
            %y(2:5:end) 
%(W*(y(2:5:end).*V))./V

            P3 = p.lambda_FU*W'*y(4:5:end) + p.lambda_FD*W*y(4:5:end) ...
               - (p.lambda_FD*(sum(W,1)'+outlet) + p.lambda_FU*sum(W,2)).*y(4:5:end);

            P4 = p.lambda_FU*W'*(y(4:5:end).*y(5:5:end))./y(4:5:end) ...
               + p.lambda_FD*W*(y(4:5:end).*y(5:5:end))./y(4:5:end) ...
               - (p.lambda_FD*(sum(W,1)'+outlet) + p.lambda_FU*sum(W,2)).*y(5:5:end);

            P4(isnan(P4)) = 0;
            dy(1:5:end) = p.alpha_M*P1./H - (p.mu_H+p.mu_W+p.gamma)*y(1:5:end);

            %dy(2:5:end) = p.rho_E * H .* y(1:5:end)./(p.omega+y(1:5:end))./V ...
            %     - (p.mu_E + p.beta_E * S./V).*y(2:5:end);

            %dy(4:5:end) = p.mu_F*(KF-y(4:5:end)) - chi.*y(4:5:end) ;

            %dy(5:5:end) = p.theta_C*S.*y(3:5:end) - (p.mu_F+chi).*y(5:5:end);

            

                 
            % dy(2:5:end) = p.rho_E * H .* y(1:5:end)./(p.omega+y(1:5:end))./V ...
            %     + p.lambda_ED * P2 - (p.mu_E + p.beta_E * S./V).*y(2:5:end);

            dy(2:5:end) = p.rho_E * H .* y(1:5:end)./(p.omega+y(1:5:end))./V ...
                - (p.mu_E + p.beta_E * S./V).*y(2:5:end);
% f
            dy(isnan(dy))=0; % because the volume is 0 in a few instances
            dy(3:5:end) = p.beta_E*(1-y(3:5:end)).*y(2:5:end) - p.mu_S*y(3:5:end);

            dy(4:5:end) = p.mu_F*(KF-y(4:5:end)) - chi.*y(4:5:end) + P3;

            dy(5:5:end) = p.theta_C*S.*y(3:5:end) + P4 - (p.mu_F+chi).*y(5:5:end);

         if any(isnan(dy))
             error();
         end
        end
    end

    
end