function y = model_diffusion(Time,par,setup,y0)

    y=odemodel(par,...
        setup.nNodes,...
        setup.WF,...#bi-directional hydrological connectivity
        setup.outlet,...
        Time,...    
        y0);
        
    %%% ODE PART
    
    function y = odemodel(par,nNodes,WF,outlet,tspan,y0)

        [~,y]=ode45(@eqs,tspan,y0);
        
        function dy=eqs(t,y)

            dy=zeros(3*nNodes,1);

            dy(3:3:end) = par.lambda_F*WF*y(3:3:end) - ...
                par.lambda_F*sum(WF,2).*y(3:3:end) - ...
                (par.mu_F+outlet).*y(3:3:end); 
         
        end
    end

    
end