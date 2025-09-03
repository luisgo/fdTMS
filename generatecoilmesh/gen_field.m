function E = gen_field(r,R,theta,phi,N,i)

    const = 1/(4*pi*10^(-7));

    temp_val = zeros(3,1);

    E = zeros(3,1);

    [P, dPx] = gen_matrix(N, cos(theta)); %good

    for L = 1:N

        for M = -L:L

            Ymll = cross_prod(M, L, r, theta, phi, P, dPx);

            scalar = (i/(2*L+1))*(r/R)^L;

            temp_val(2,1) = Ymll(2)*scalar + temp_val(2,1);

            temp_val(3,1) = Ymll(3)*scalar + temp_val(3,1);

        end

    end

    E(2,1) = const*temp_val(2,1);

    E(3,1) = const*temp_val(3,1);

end



function Ymll = cross_prod(M, L, r, theta, phi, P, dPx)

    dYml = gen_gradient(M,L,theta,phi, P,dPx);

    rvector = zeros(3,1);

    rvector(1,1) = (sqrt(L*(L+1)))^(-1) * r;

    Ymll = cross(rvector, dYml);

end



function dYml = gen_gradient(M,L,theta,phi,P,dPx)

    dYml = zeros(3,1);

    if M > 0

        const = sqrt((2*L+1)*factorial(L-M)/(2*pi*factorial(L+M)));

        dYml(2,1) = const*dPx(L,M+1)*(-sin(theta))*cos(M*phi);

        dYml(3,1) = const*P(L,M+1)*(-M*sin(M*phi));

    elseif M == 0

        dYml(2,1) = sqrt((2*L+1)/(4*pi))*dPx(1,1)*(-sin(theta));

    else

        const = -sqrt((2*L+1)*factorial(L+M)/(2*pi*factorial(L-M)));

        dYml(2,1) = const*dPx(L,-M+1)*(-sin(theta))*sin(M*phi);

        dYml(3,1) = const*P(L,-M+1)*(M*cos(M*phi));

    end

end



function [P, dPx] = gen_matrix(N, x)

    P = zeros(N, N+1);

    dPx = zeros(N, N+1);

    for M = 0:N

        if M == 0

            base0 = 1;

            base1 = x;

            P(1,M+1) = base0;

            P(2,M+1) = base1;

        elseif M == 1

            base0 = -sqrt(1-x^2);

            base1 = 3*x*base0;

        else 

            base0 = (-1)^M*prod(1:2:(2*M-1))*sqrt((1-x^2)^M);

            base1 = (2*M+1)*x*base0;

        end

        for L = 1:N

            if L == M 

                P(L,M+1) = base0;

            elseif L == M+1 && M~=0

                P(L,M+1) = base1;

            elseif L >= 3 

                P(L, M+1) = (x*(2*L-1)*P(L-1,M+1) - (L+M-1)*P(L-2))/(L-M);

            end 

        end    

    end

    for M = 0:N

        if M == 0

            base0 = 0;

            base1 = 1;

            P(1,M+1) = base0;

            P(2,M+1) = base1;

        elseif M == 1

            base0 = x*(1-x^2)^(-1/2);

            base1 = (2*M+1)*P(M,M+1) + base0;

        else 

            base0 = (-1)^M*prod(1:2:(2*M-1))*(x*M*(1-x^2)^(M/2+1));

            base1 = (2*M+1)*P(M,M+1) + base0;

        end

        for L = 1:N

            if L == M 

                dPx(L,M+1) = base0;

            elseif L == M+1 && M~=0

                dPx(L,M+1) = base1;

            elseif L >= 3 

                dPx(L,M+1) = ((2*L-1)*P(L-1,M+1) + x*(2*L-1)*dPx(L-1,M+1) + (L+M-1)*dPx(L-2,M+1))/(L-M);

            end 

        end

    end             

end

