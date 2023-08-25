function y = model(Time,par,setup,y0)

    y=odemodel( ...
        par,...
        setup.nNodes,... #number of nodes
        setup.N_F,...
        setup.Cf,...#caught fish quota
        setup.M,... #fish market matrix
        setup.W,...#bi-directional hydrological connectivity
        setup.outlet,...
        Time,...    
        y0 ...
        );
        
    %%% ODE PART
    
    function y = odemodel(par,nNodes,N_F,Cf,M,W,outlet,tspan,y0)

        [~,y]=ode45(@eqs,tspan,y0);
        
        function dy=eqs(t,y)

            dy=zeros(1*nNodes,1);

            dy(3:3:end) = par.mu_F*(N_F-y) +...
            (par.lambda_FU*W+par.lambda_FD*W')*(y)./N_F - ...
            (par.lambda_FD*sum(W,2)+par.lambda_FU*sum(W,1)').*y(3:3:end);
            
            
            
            %- ...
               % par.lambda_F*sum(WF,2).*y(3:3:end) - ...
                %(par.mu_F+outlet+Cf).*y(3:3:end); 
         
        end
    end

    
end