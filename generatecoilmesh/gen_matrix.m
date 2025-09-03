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

                P(L, M+1) = (x*(2*L-1)*P(L-1,M+1) - (L+M-1)*P(L-2,M+1))/(L-M);

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