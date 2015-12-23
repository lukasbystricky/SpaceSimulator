function [P, Pd] = associated_legendre(N, x, norm_type)

%% Evaluates all fully normalized associated Legendre polynomials 
% of degree and order less than or equal to N-1, and their first derivatives at x.
%
% For n greater than or equal to m:
% P(m+1,n+1) = P_n^m(x)
% Pd(m+1, n+1) = d/dx P_n^m(x)

P   = zeros(N, N);
Pd = zeros(N, N);
full_derivs = zeros(N+1,N);

%% compute P^n_0
[L, Ld] = legendre_order0(N,x);
P(1,:) = L;
Pd(1,:) = Ld;

%% calculate all derivatives up to order N of of P^n_0
full_derivs(1,:) = Pd(1,:);
for j = 2:N+1
    for i = 3:N
        n = i - 2;
        full_derivs(j,i) = ((2*n + 1)*(j*full_derivs(j-1,i-1)+ x*full_derivs(j,i-1)) ...
			- n*full_derivs(j,i-2))/(n + 1);
    end
end

%% compute P^m_n
for m = 2:N
    for n = m:N
        if ~(strcmp(norm_type, 'positive') || strcmp(norm_type, 'egm96'))
            factor = (-1)^(m-1);
        else
            factor = 1;
        end
        P(m,n) = factor*(1 - x^2)^((m-1)/2)*full_derivs(m-1,n);
        Pd(m,n) = -factor*((m-1)/2)*(1 - x^2)^(((m-1)/2)-1)*2*x + ...
            factor*(1-x^2)^((m-1)/2)*full_derivs(m,n);
    end
end

%% normalize values
if ~strcmp(norm_type, 'none')
    for n = 0:N-1
        for m = 0:N-1
            if (n >= m)
                switch norm_type
                    
                    case 'standard'
                        factor = ((-1)^m)*sqrt((n+0.5)*factorial(n-m)/factorial(n+m));
                        
                    case 'power'
                        factor = sqrt(factorial(n-m)/factorial(n+m));
                        
                    case 'positive'
                        factor = 1;
                      
                    case 'egm96'
                        factor = (2*n + 1)*factorial(n - m)/factorial(n + m);
                        
                        if (m ~= 0)
                            factor = factor*2;
                        end
                        
                        factor = sqrt(factor);
                        
                    otherwise
                        disp('invalid normalization type, valid inputs are: none, standard and power');
                        break;
                        
                end
                P(m+1, n+1) = P(m+1, n+1)*factor;
                Pd(m+1, n+1) = Pd(m+1, n+1)*factor;
            end
        end
    end
end
end


function [L, Ld] = legendre_order0(N,x)

L = zeros(1,N);
Ld = zeros(1,N);

L(1) = 1;
L(2) = x;
Ld(1) = 0;
Ld(2) = 1;

for i = 3:N
    n = i - 2;
    L(i) = ((2*n + 1)*x*L(i-1) - n*L(i-2))/(n+1);
    Ld(i) = ((2*n + 1)*(L(1,i-1) + x*Ld(1,i-1)) - n*Ld(1,i-2))/(n+1);
end
end
