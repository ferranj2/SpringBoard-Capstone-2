function [TtT,PtP,MFP] = MASSFP_NR(Tt,f,M)
C0 = (f*+7.3816638e-2 +2.5020051e-1)/(1+f);
C1 = (f*+1.2258630e-3 -5.1536879e-5)/(1+f); D1 = C1/2;
C2 = (f*-1.3771901e-6 +6.5519486e-8)/(1+f); D2 = C2/3;
C3 = (f*+9.9686793e-10 -6.7178376e-12)/(1+f); D3 = C3/4; T3 = C3*2;
C4 = (f*-4.2051104e-13 -1.5128259e-14)/(1+f); D4 = C4/5; T4 = C4*3;
C5 = (f*+1.0212913e-16 +7.6215767e-18)/(1+f); D5 = C5/6; T5 = C5*4;
C6 = (f*-1.3335668e-20 -1.4526770e-21)/(1+f); D6 = C6/7; T6 = C6*5;
C7 = (f*+7.2678710e-25 +1.0115540e-25)/(1+f); D7 = C7/8; T7 = C7*6;
h_ref = (f*30.58153 -1.7558886)/(1+f); 
[~,ht,Prt,~,~,R,~,~] = FAIR(1,Tt,f,'BE');
const1 = (M^2)*R/2;
const2 = ((M*R)^2)/2;
k = 0;
k_max = 100;
diff = 0.5;
if M > 1
    T = Tt*0.7;
else
    T = Tt*0.95;
end
while diff > 0.000001 && k < k_max
    Cp = C0 + C1*T + C2*(T^2) + C3*(T^3) + C4*(T^4) + C5*(T^5) + C6*(T^6) + C7*(T^7);
    h = h_ref + C0*T + D1*(T^2) + D2*(T^3) + D3*(T^4) + D4*(T^5) + D5*(T^6) + D6*(T^7) + D7*(T^8);
    h_iter = h - ht + const1*Cp*T/(Cp-R);
    h_prime = Cp + const2*(power(Cp,2)/R - 2*Cp + C0 - C2*(T^2) - T3*(T^3) - T4*(T^4) - T5*(T^5) - T6*(T^6) - T7*(T^7))/power(Cp-R,2);
    Tnew = T - h_iter/h_prime;
    k = k + 1;
    diff = abs(Tnew-T);
    T = Tnew;
end
TtT = Tt/T;
[~,~,Pr,~,~,~,gamma,~] = FAIR(1,T,f,'BE');
PtP = Prt/Pr;
MFP = (M/PtP)*sqrt(gamma*TtT*32.174/R);
%error = ht - h - ((M*a)^2)/(2*32.174*778.16)
end