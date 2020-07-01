%TO-DO
%- Collect all thermo outputs.

%--------------------------------------------------------------------
%Title: "Mattingly's LBTF Parametric Cycle Analysis w/ VSH in MATLAB"
%By: Jesus Ferrand
%Credits: Jack D. Mattingly
%--------------------------------------------------------------------
%INFO:
function R = ONX_LBTF_VSH(I)
R = struct('A',I.A,'A_prime',0,'A4_A45',0,'A45_A6',0,'A16_A6',0,'A8_A6_dry',0,...
    'A6_A6A',0,'A9_A8_off',0,'A9_A8_on',0,'f',0,'f_AB',0,'f_6A',0,'f_0_on',0,'f_0_off',0,'F_m_dot_on',0,...
    'F_m_dot_off',0,'MFP_4',0,'MFP_45',0,'MFP_6',0,'MFP_6A',0,...
    'pi_ABdry',0,'pi_AB',I.pi_AB,'PtP_9_off',0,'PtP_9_on',0,...
    'pi_c',I.pi_c,...
    'pi_cL',I.pi_cL,'pi_cH',0,'pi_d',0,'pi_f',I.pi_f,'pi_M',0,...
    'pi_M_ideal',0,'pi_r',0,'pi_tH',0,'pi_tL',0,'eta_cH',0,'eta_cL',0,...
    'eta_tH',0,'eta_tL',0,'eta_f',0,'eta_P_on',0,'eta_P_off',0,...
    'eta_TH_on',0,'eta_TH_off',0,'m_dot',I.m_dot,'M_4',1,'M_45',1,...
    'M_6',I.M_6,'M_16',0,'M_6A',0,'M_8_off',0,'M_8_on',0,'M_9_on',0,'M_9_off',0,...
    'P_0',0,...
    'h_0',0,...
    'Pt16_Pt6',0,...
    'S_on',0,...
    'Cp_6A',0,'gamma_6A',0,...
    'S_off',0,'tau_cH',0,'tau_cL',0,'tau_f',0,'tau_m1',0,'tau_m2',0,...
    'tau_M',0,'tau_tH',0,'tau_tL',0,'tau_lambda',0,'tau_r',0,'theta_0',0,...
    'T9_T0_off',0,'T9_T0_on',0,'V9_V0_off',0,'V9_V0_on',0,'M9_M0_off',0,'M9_M0_on',0,...
    'P_TOL',0,'P_TOH',0,'T_0',0,'Tt_45',0,'Tt_4',I.Tt_4,'Tt_7',I.Tt_7,'e1',I.e1,'e2',I.e2,'B',I.B,...
    'lmbme1me2',1-I.B-I.e1-I.e2);
g_c = 32.174;
R.pi_ABdry = 1-(1-I.pi_AB)/2;
[R.T_0,~,R.P_0,rho_0] = ATMOSphere(I.Alt,'BE',I.day);
[~,R.h_0,Pr_0,Phi_0,Cp_0,R_0,gamma_0,a_0] = FAIR(1,R.T_0,0,'BE');
R.P_TOL = I.CTOL*I.m_dot*R.h_0; %Low-Pressure Power Off-Take [Btu/s]
R.P_TOH = I.CTOH*I.m_dot*R.h_0; %High-Pressure Power Off-Take [Btu/s]
V_0 = I.M_0*a_0; %ONX Ambient Velocity. [ft/s]
ht_0 = R.h_0 + (V_0^2)/(2*g_c*778.16); %Stagnation Enthalpy at Ambient.
R.tau_r = ht_0/R.h_0; %Ram Heating.
R.theta_0 = R.tau_r*R.T_0/518.67; %ONX Throttle Ratio.
[Tt_0,ht_0,Prt_0,Phit_0,~,~,~,~] = FAIR(2,ht_0,0,'BE');
R.pi_r = Prt_0/Pr_0; %Ram Pressure Ratio .
eta_R_spec = MIL_E_5008B(I.M_0);
R.pi_d = I.pi_d_max*eta_R_spec; %Diffuser Pressure Ratio.
ht_2 = ht_0; %Stagnation Enthalpy entering Fan. (Adiabatic Assumption)
Prt_2 = Prt_0; %Reduced Pressure entering Fan. (Adiabatic Assumption)
Prt_13 = Prt_2*power(I.pi_f,1/I.e_f);
[Tt_13,ht_13,~,Phit_13,~,~,~,~] = FAIR(3,Prt_13,0,'BE');
R.tau_f = ht_13/ht_2; %Fan Enthalpy ratio.
Prt_13i = Prt_2*I.pi_f;
[~,ht_13i,~,~,~,~,~,~] = FAIR(3,Prt_13i,0,'BE');
R.eta_f = (ht_13i-ht_2)/(ht_13-ht_2);
Prt_25 = Prt_2*power(I.pi_cL,1/I.e_cL);
[Tt_25,ht_25,~,Phit_25,~,~,~,~] = FAIR(3,Prt_25,0,'BE');
R.tau_cL = ht_25/ht_2;
Prt_25i = Prt_2*I.pi_cL;
[~,ht_25i,~,~,~,~,~,~] = FAIR(3,Prt_25i,0,'BE');
R.eta_cL = (ht_25i-ht_2)/(ht_25-ht_2);
R.pi_cH = I.pi_c/I.pi_cL;
Prt_3 = Prt_25*power(R.pi_cH,1/I.e_cH);
[Tt_3,ht_3,~,Phit_3,~,~,~,~] = FAIR(3,Prt_3,0,'BE');
R.tau_cH = ht_3/ht_25;
Prt_3i = Prt_25*R.pi_cH;
[~,ht_3i,~,~,~,~,~,~] = FAIR(3,Prt_3i,0,'BE');
R.eta_cH = (ht_3i-ht_25)/(ht_3-ht_25);
R.f = FUEL_NR(ht_3,I.eta_b,I.h_PR,I.Tt_4);
lmbmelme2tlpf = (1-R.B-R.e1-R.e2)*(1+R.f); %Optimization Step.
[~,ht_4,Prt_4,Phit_4,~,~,~,~] = FAIR(1,R.Tt_4,R.f,'BE');
[Tt4_T4,Prt4_Pr4,R.MFP_4] = MASSFP_NR(R.Tt_4,R.f,R.M_4);
R.tau_lambda = ht_4/R.h_0;
R.tau_m1 = (lmbmelme2tlpf+I.e1*R.tau_r*R.tau_cL*R.tau_cH/R.tau_lambda)/(lmbmelme2tlpf+I.e1);
R.tau_tH = 1 - (R.tau_r*R.tau_cL*(R.tau_cH-1)+(1+I.A)*I.CTOH/I.eta_mPH)/((I.eta_mH*R.tau_lambda)*lmbmelme2tlpf+I.e1*R.tau_r*R.tau_cL*R.tau_cH/R.tau_lambda);
ht_41 = ht_4*R.tau_m1;
f_41 = R.f/(1+R.f+I.e1/R.lmbme1me2);
[Tt_41,~,Prt_41,Phit_41,~,~,~,~] = FAIR(2,ht_41,f_41,'BE');
ht_44 = ht_41*R.tau_tH;
[Tt_44,~,Prt_44,Phit_44,~,~,~,~] = FAIR(2,ht_44,f_41,'BE');
R.pi_tH = power(Prt_44/Prt_41,1/I.e_tH);
Prt_44i = R.pi_tH*Prt_41;
[~,ht_44i,~,~,~,~,~,~] = FAIR(3,Prt_44i,f_41,'BE');
R.eta_tH = (ht_41-ht_44)/(ht_41-ht_44i);
R.tau_m2 = (lmbmelme2tlpf+I.e1+I.e2*(R.tau_r*R.tau_cL*R.tau_cH/(R.tau_lambda*R.tau_m1*R.tau_tH)))/(lmbmelme2tlpf+I.e1+I.e2);
ht_45 = ht_44*R.tau_m2;
f_45 = R.f/(1+R.f+(I.e1+I.e2)/R.lmbme1me2);

[R.Tt_45,~,Prt_45,Phit_45,~,~,~,~] = FAIR(2,ht_45,f_45,'BE');
[TtT45,PtP45,R.MFP_45] = MASSFP_NR(R.Tt_45,f_45,R.M_45);
R.A4_A45 = R.pi_tH*sqrt(R.Tt_4/R.Tt_45)*(R.MFP_45/R.MFP_4)/(1+(I.e1+I.e2)/lmbmelme2tlpf);
R.tau_tL = 1 - (R.tau_r*((R.tau_cL-1)+I.A*(R.tau_f-1))+(1+I.A)*I.CTOL/I.eta_mPL)/  (I.eta_mL*R.tau_lambda*R.tau_tH*(lmbmelme2tlpf+(I.e1+I.e2/R.tau_tH)*(R.tau_r*R.tau_cL*R.tau_cH/R.tau_lambda)));
ht_5 =ht_45*R.tau_tL;
[Tt_5,~,Prt_5,Phit_5,~,~,~,~] = FAIR(2,ht_5,f_45,'BE');
R.pi_tL = power(Prt_5/Prt_45,1/I.e_tL);
Prt_5i = R.pi_tL*Prt_45;
[~,ht_5i,~,~,~,~,~,~] = FAIR(3,Prt_5i,f_45,'BE');
R.eta_tL = (ht_45-ht_5)/(ht_45-ht_5i);

%Mixer Calculations
R.A_prime = I.A/(lmbmelme2tlpf+I.e1+I.e2);
f_6 = f_45;
R.f_6A = f_6/(1+R.A_prime);
ht_16 = ht_13;
Prt_16 = Prt_13;
ht_6 = ht_5; 
ht_6A = (ht_6+R.A_prime *ht_16)/(1+R.A_prime);
R.tau_M = ht_6A/ht_6;
R.Pt16_Pt6 = I.pi_f/(I.pi_c*I.pi_b*R.pi_tH*R.pi_tL);
Tt_6 = Tt_5;

Tt_16 = Tt_13;
[Tt6_T6,Pt6_P6,R.MFP_6] = MASSFP_NR(Tt_6,f_45,I.M_6);
R.A45_A6 = R.pi_tL*sqrt(R.Tt_45/Tt_6)*(R.MFP_6/R.MFP_45);
PtP16 = Pt6_P6*R.Pt16_Pt6;
Pr_16 = Prt_16/PtP16;
[T_16,h_16,~,Phi_16,Cp_16,R_16,gamma_16,a_16] = FAIR(3,Pr_16,0,'BE');
V_16 = sqrt(2*g_c*(ht_16-h_16)*778.16);
R.M_16 = V_16/a_16;
[TtT16,PtP16,MFP_16] = MASSFP_NR(Tt_16,0,R.M_16);
R.A16_A6 = R.A_prime*sqrt(Tt_16/Tt_6)*(R.MFP_6/MFP_16)/R.Pt16_Pt6;
R.A6_A6A = 1/(1+R.A16_A6);
T_6 = Tt_6/Tt6_T6;
[~,h_6,Pr_6,Phi_6,Cp_6,R_6,gamma_6,a_6] = FAIR(1,T_6,f_6,'BE');
[R.M_6A,R.MFP_6A,Tt_6A,R.gamma_6A,R.Cp_6A] = MIXER_VSH(I.M_6,R.M_16,gamma_16,gamma_6,R_6,R.A_prime,R.A16_A6,T_6,R.f_6A,ht_6A);
R.pi_M_ideal = (1+R.A_prime)*sqrt(R.tau_M)*R.A6_A6A*R.MFP_6/R.MFP_6A;
R.pi_M = R.pi_M_ideal*I.pi_M_max;
R.f_AB = FUEL_NR(ht_6A,I.eta_AB,I.h_PR,I.Tt_7);

%Afterburner OFF
R.f_0_off = R.lmbme1me2*R.f/(1+I.A); % Overall fuel-air ratio (AB = off). [Dimensionless]
[~,ht_7_off,Prt_7_off,Phit_7__off,~,~,~,~] = FAIR(1,Tt_6A,R.f_0_off,'BE');
ht_9_off = ht_7_off;
Prt_9_off = Prt_7_off;
R.PtP_9_off = I.P0_P9*R.pi_r*R.pi_d*R.pi_cL*R.pi_cH*I.pi_b*R.pi_tH*R.pi_tL*R.pi_M*R.pi_ABdry*I.pi_n;
Pr_9_off = Prt_9_off/R.PtP_9_off;
[T_9_off,h_9_off,~,Phi_9_off,Cp_9_off,R_9_off,gamma_9_off,a_9_off] = FAIR(3,Pr_9_off,R.f_0_off,'BE');
V_9_off = sqrt(2*g_c*778.16*(ht_9_off-h_9_off));
R.M_9_off = V_9_off/a_9_off;
if R.M_9_off > 1
    R.M_8_off = 1;
else
    R.M_8_off = R.M_9_off;
end
[TtT_9_off,~,MFP_9_off] = MASSFP_NR(Tt_6A,R.f_6A,R.M_9_off);
[TtT_8_off,PtP_8_off,MFP_8_off] = MASSFP_NR(Tt_6A,R.f_6A,R.M_8_off);
R.A8_A6_dry = R.MFP_6*(1+R.A_prime)*sqrt(Tt_6A/Tt_6)/(R.pi_M*R.pi_ABdry*MFP_8_off);
R.A9_A8_off = MFP_8_off/(MFP_9_off*I.pi_n);
R.T9_T0_off =  T_9_off/R.T_0;
R.V9_V0_off = V_9_off/V_0;
R.M9_M0_off = R.M_9_off/I.M_0;
R.F_m_dot_off = (a_0/g_c)*((1+R.f_0_off - I.B/(1+I.A))*(V_9_off/a_0 + R_9_off*T_9_off*a_0*(1-I.P0_P9)/(R_0*R.T_0*V_9_off*gamma_0)) - I.M_0); % Uninstalled Specific Thrust (AB = on). [lbf/lbm*s]
R.S_off = 3600*R.f_0_off/R.F_m_dot_off; %Uninstalled TSFC (AB = on). [1/hr]
R.eta_P_off = (2*g_c*I.M_0*R.F_m_dot_off/a_0)/((1 + R.f_0_off - I.B/(1+I.A))*((V_9_off/a_0)^2)-I.M_0^2);
R.eta_TH_off =(((1+R.f_0_off-I.B/(1+I.A))*(V_9_off^2) - V_0^2)/(2*g_c*778.16)+(I.CTOL+I.CTOH)*R.h_0)/(R.f_0_off*I.h_PR);

%Afterburner ON
R.f_0_on = (R.lmbme1me2*R.f + (R.f_AB*(1+I.A-I.B)))/(1+I.A); % Overall fuel/air ratio (AB = on). [Dimensionless]
[~,ht_7,Prt_7,Phit_7,~,~,~,~] = FAIR(1,I.Tt_7,R.f_0_on,'BE');
ht_9 = ht_7;
Prt_9 = Prt_7;
R.PtP_9_on = I.P0_P9*R.pi_r*R.pi_d*R.pi_cL*R.pi_cH*I.pi_b*R.pi_tH*R.pi_tL*R.pi_M*I.pi_AB*I.pi_n;
Pr_9 = Prt_9/R.PtP_9_on;
[T_9,h_9,~,Phi_9,Cp_9,R_9,gamma_9,a_9] = FAIR(3,Pr_9,R.f_0_on,'BE');
V_9 = sqrt(2*g_c*778.16*(ht_9-h_9));
R.M_9_on = V_9/a_9;
if R.M_9_on > 1
    R.M_8_on = 1;
else
    R.M_8_on = R.M_9_on;
end
[TtT_9_on,~,MFP_9_on] = MASSFP_NR(I.Tt_7,R.f_0_on,R.M_9_on);
[TtT_8_on,PtP_8_on,MFP_8_on] = MASSFP_NR(I.Tt_7,R.f_0_on,R.M_8_on);
R.A9_A8_on = MFP_8_on/(MFP_9_on*I.pi_n);
R.T9_T0_on =  T_9/R.T_0;
R.V9_V0_on = V_9/V_0;
R.M9_M0_on = R.M_9_on/I.M_0;
R.F_m_dot_on = (a_0/g_c)*((1+R.f_0_on - I.B/(1+I.A))*(V_9/a_0 + R_9*T_9*a_0*(1-I.P0_P9)/(R_0*R.T_0*V_9*gamma_0)) - I.M_0); % Uninstalled Specific Thrust (AB = on). [lbf/lbm*s]
R.S_on = 3600*R.f_0_on/R.F_m_dot_on; %Uninstalled TSFC (AB = on). [1/hr]
R.eta_P_on = (2*g_c*I.M_0*R.F_m_dot_on/a_0)/((1 + R.f_0_on - I.B/(1+I.A))*((V_9/a_0)^2)-I.M_0^2);
R.eta_TH_on =(((1+R.f_0_on-I.B/(1+I.A))*(V_9^2) - V_0^2)/(2*g_c*778.16)+(I.CTOL+I.CTOH)*R.h_0)/(R.f_0_on*I.h_PR);
end