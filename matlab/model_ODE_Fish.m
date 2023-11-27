function y = model_ODE_Fish(Time,F,H)

    W = zeros(size(H));
    startNode = 1; % Choose any starting node
    visited = false(size(H,1), 1);
    W = calculateWRecursive(startNode, H, F, W, visited,1)



    [~,y]=ode45(@eqs,Time,F);
    
    function dy=eqs(t,y)

        dy = W*y - sum(W,1)'.*y;
            
    end
    close all
    figure(); plot(y)
    
end