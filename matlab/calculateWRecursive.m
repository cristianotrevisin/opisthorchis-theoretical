function W = calculateWRecursive(node, H, F, W, visited,val)
    
    connectedNodes = find(H(node, :) == 1);
    a = find(W(connectedNodes,node)>0);
    
    % Assign ones to first node
    W(connectedNodes,node) = val;

    visited(node) = true;
    % Start recursive computation (exploration along each branch)
    for i = connectedNodes
        if ~visited(i)
            % Calculate W(i, node) and W(node, i)
            W(node, i) = W(i,node) * F(node) / F(i);
            val = W(i,node) * F(node) / F(i);

            W = calculateWRecursive(i, H, F, W, visited,val);
            visited(i) = true;
        end
    end
end