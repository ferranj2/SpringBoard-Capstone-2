function [pi_tH,tau_tH,Tt_45] = TURBC(Tt_4,ht_4,f,A4_A45,M_4,M_45,eta_tH,Tt_45R,ht_3,lmbme1me2,e1,e2)
[TtT4,PtP4,MFP_4] = MASSFP_NR(Tt_4,f,M_4);
m_f = f*lmbme1me2;
m_4 = (1+f)*lmbme1me2;
m_41 = m_4+e1;
m_45 = m_41+e2;
f_41 = m_f/(m_41-m_f);
f_45 = m_f/(m_45-m_f);
ht_41 = (m_4*ht_4+e1*ht_3)/(m_4+e1);
[Tt_41,~,Prt_41,Phit_41,~,~,~,~] = FAIR(2,ht_41,f_41,'BE');
Tt_45 = Tt_45R;
k = 0;
k_max = 100;
diff = 0.5;
while diff > 0.01 && k < k_max
[TtT45,PtP45,MFP_45] = MASSFP_NR(Tt_45,f_45,M_45);
pi_tH = sqrt(Tt_45/Tt_4)*A4_A45*m_45*MFP_4/(m_4*MFP_45);
%[~,h,Pr,Phi,~,~,~,~] = FAIR(1,Tt_45R,f_45,'BE'); %??????
Prt_44i = pi_tH*Prt_41;
[~,ht_44i,Pr_44i,Phi_44i,~,~,~,~] = FAIR(3,Prt_44i,f_41,'BE');
ht_44 = ht_41-eta_tH*(ht_41-ht_44i);
tau_tH = ht_44/ht_41;
ht_45 = (m_41*ht_44+e2*ht_3)/(m_41+e2);
[Tt_45n,h,Pr,Phi,~,~,~,~] = FAIR(2,ht_45,f_45,'BE');
diff = abs(Tt_45n-Tt_45);
k = k+1;
Tt_45 = Tt_45n;
end
end