function f = FUEL_NR(enthalpy,eta_b,h_PR,Tt)
A0 = +2.5020051e-1;B0 = +7.3816638e-2;
A1 = -5.1536879e-5;B1 = +1.2258630e-3;
A2 = +6.5519486e-8;B2 = -1.3771901e-6;
A3 = -6.7178376e-12;B3 = +9.9686793e-10;
A4 = -1.5128259e-14;B4 = -4.2051104e-13;
A5 = +7.6215767e-18;B5 = +1.0212913e-16;
A6 = -1.4526770e-21;B6 = -1.3335668e-20;
A7 = +1.0115540e-25;B7 = +7.2678710e-25;
Ah_ref = -1.7558886;Bh_ref = 30.58153; % [Btu/lbm]

k = 0;
k_max = 100;
error = 0.5;

Asum = sum([Ah_ref,A0*Tt,A1*(Tt^2)/2,A2*(Tt^3)/3,A3*(Tt^4)/4,A4*(Tt^5)/5,A5*(Tt^6)/6,A6*(Tt^7)/7,A7*(Tt^8)/8]);
Bsum = sum([Bh_ref,B0*Tt,B1*(Tt^2)/2,B2*(Tt^3)/3,B3*(Tt^4)/4,B4*(Tt^5)/5,B5*(Tt^6)/6,B6*(Tt^7)/7,B7*(Tt^8)/8]);
f = 0.031;% Initial Guess
while error > 0.0000000000000001 && k <= k_max
    num = Asum - enthalpy + f*(Bsum-enthalpy);
    den = eta_b*h_PR-Asum+f*(eta_b*h_PR-Bsum);
    f_iter = num/den - f;
    f_prime = ((Bsum-enthalpy)*den-num*(eta_b*h_PR-Bsum))/(den^2) - 1;
    f_new = f - f_iter/f_prime;
    error = abs(f_new-f);
    k = k+1;
    f = f_new;
end
end