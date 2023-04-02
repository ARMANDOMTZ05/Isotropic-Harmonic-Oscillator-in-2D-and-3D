% Codigo para evaluar los polinmios de Laguerre

function Poli = Laguerre_pol(x,a,n)
    poli_1 = zeros(1,length(x));
    poli_2 = zeros(1,length(x));
    poli_3 = zeros(1,length(x));

    if n == 0
        poli_1 = 1;
    elseif n == 1
        poli_2 = 1+a-x;
    else
        poli_3 = ((2*n+a-1-x)/n).*Laguerre_pol(x,a,n-1) - ((n+a-1)/n).*Laguerre_pol(x,a,n-2);
    end
    
    Poli = poli_1 + poli_2 + poli_3;
end