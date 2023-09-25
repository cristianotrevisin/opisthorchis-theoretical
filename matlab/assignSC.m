function M = assignSC(MSC,X)
    M = zeros(size(MSC));
    for iz = 1:length(X)
        M(MSC==iz) = X(iz);
    end
end