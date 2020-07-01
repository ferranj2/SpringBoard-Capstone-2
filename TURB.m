function [pi_tL,tau_tL,Tt_e] = TURB(Tt_i,f,Ai_Ae,M_i,M_e,eta_tL,Tt_eR)
[~,ht_i,Prt_i,Phit_i,~,~,~,~] = FAIR(1,Tt_i,f,'BE');
[TtTi,PtPi,MFP_i] = MASSFP_NR(Tt_i,f,M_i);
Tt_e = Tt_eR;
k = 0;
k_max = 100;
diff = 0.5;
while diff > 0.01 && k < k_max
    [TtTe,PtPe,MFP_e] = MASSFP_NR(Tt_e,f,M_e);
    pi_tL = MFP_i*Ai_Ae*sqrt(Tt_e/Tt_i)/MFP_e;
    Prt_ei = pi_tL*Prt_i;
    [~,ht_ei,Pr_ei,Phi_ei,~,~,~,~] = FAIR(3,Prt_ei,f,'BE');
    ht_e = ht_i - eta_tL*(ht_i-ht_ei);
    tau_tL = ht_e/ht_i;
    [Tt_en,~,~,~,~,~,~,~] = FAIR(2,ht_e,f,'BE');
    diff = abs(Tt_en-Tt_e);
    Tt_e = Tt_en;
    k = k+1;
end
end