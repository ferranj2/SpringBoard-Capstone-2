function [M_6A,MFP_6A,Tt_6A,gamma_6A,Cp_6A] = MIXER_VSH(M_6,M_16,gamma_16,gamma_6,R_6,alpha_prime,A16_A6,T_6,f_6A,ht_6A)
C = sqrt(R_6*T_6/gamma_6)*((1+gamma_6*M_6^2)+A16_A6*(1+gamma_16*M_16^2))/(M_6*(1+alpha_prime));
M_p = M_6; % Initially guess that M_6A = M_6.
error = 1;
k = 0;
[Tt_6A,~,Prt_6A,~,~,~,~,~] = FAIR(2,ht_6A,f_6A,'BE');
[TtT6A,~,MFP_6A] = MASSFP_NR(Tt_6A,f_6A,M_p);
while abs(error) > 0.00000001 && k <= 100
    T_6A = Tt_6A/TtT6A;
    [~,h_6A,Pr_6A,Phi_6A,Cp_6A,R_6A,gamma_6A,a_6A] = FAIR(1,T_6A,f_6A,'BE');
    f = sqrt(R_6A*T_6A/gamma_6A)*(1+gamma_6A*M_p^2)/M_p - C; % Non-Linear function to be iterated.
    d = (-1/(M_p^2) + gamma_6A*M_p)*sqrt(R_6A*T_6A/gamma_6A); % Derivative of the Non-Linear function.
    M_6A = M_p - f/d; % New approximation resulting from the iteration.
    M_p = M_6A;
   [TtT6A,~,MFP_6A] = MASSFP_NR(Tt_6A,f_6A,M_p);
    error = abs(M_6A-M_p); % Update the error measured.
    k = k+1;
end
end