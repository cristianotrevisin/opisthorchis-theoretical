function P = zipf(rank, expn, minP)

    H = sum(1./(rank.^expn));
    
    
    P = 1./rank.^expn./H;
    P = P'/min(P)*minP;
    
end