function P2 = chebyshev_polynomial(x, n)

P0 = 1;
P1 = x;

if n == 0
    P2 = P0;
else if n == 1
        P2 = P1;
    else
        
        for i = 2:n
            P2 = 2*x*P1 - P0;
            P1 = P2;
            P0 = P1;
        end
    end
end
    
