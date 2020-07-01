function [TtT,PtP,M] = INV_MASSP_NR(mode,Tt,f,MFP)
if mode == 4 %Supersonic solution Desired
    T = Tt*0.6;
else %Subsonic solution Desired.
    T = Tt*0.95;
end
k = 0;
k_max = 100;
diff = 0.5;
%Pre-computed values:
[~,ht,Prt,~,~,R,~,~] = FAIR(1,Tt,f,'BE');
C0 = (f*+7.3816638e-2 +2.5020051e-1)/(1+f); 
C1 = (f*+1.2258630e-3 -5.1536879e-5)/(1+f); D1 = C1/2;
C2 = (f*-1.3771901e-6 +6.5519486e-8)/(1+f); D2 = C2/3; S2 = C2/2; 
C3 = (f*+9.9686793e-10 -6.7178376e-12)/(1+f); D3 = C3/4; S3 = C3/3;  
C4 = (f*-4.2051104e-13 -1.5128259e-14)/(1+f); D4 = C4/5; S4 = C4/4; 
C5 = (f*+1.0212913e-16 +7.6215767e-18)/(1+f); D5 = C5/6; S5 = C5/5; 
C6 = (f*-1.3335668e-20 -1.4526770e-21)/(1+f); D6 = C6/7; S6 = C6/6; 
C7 = (f*+7.2678710e-25 +1.0115540e-25)/(1+f); D7 = C7/8; S7 = C7/7; 
h_R = (f*30.58153 -1.7558886)/(1+f);
Phi_R = (f*0.6483398+0.0454323)/(1+f);
Phi_R0 =  Phi_R - 1.5784416522122042;
C = R*(MFP*Prt)/sqrt(64.348*Tt);
while diff > 0.000001 && k < k_max    
    Cp = C7*T^7 + C6*T^6 + C5*T^5 + C4*T^4 + C3*T^3 + C2*T^2 + C1*T + C0;
    Pr = exp((Phi_R0 + C0*log(T) + S2*T^2 + S3*T^3 + S4*T^4 + S5*T^5 + S6*T^6 + S7*T^7 + C1*T)/R);
    htmh = -(D7*T^8 + D6*T^7 + D5*T^6 + D4*T^5 + D3*T^4 + D2*T^3 + D1*T^2 + C0*T + h_R - ht);
    dhdT2 = (htmh/T^2)^(1/2);
    f = Pr*dhdT2 - C;
    f_der = Pr*((dhdT2*(C0/T + C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4 + C6*T^5 + C7*T^6))/R - ((Cp/T^2 + (2*htmh)/T^3))/(2*dhdT2));
    Tn = T-f/f_der;
    diff = abs(Tn-T);
    T = Tn;
    k = k+1;
end
TtT = Tt/T;
PtP = Prt/Pr;
M = sqrt(2*(htmh)/(R*T*Cp/(Cp-R)));
end