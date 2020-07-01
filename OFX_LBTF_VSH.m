clc;
%clear;
%close all;
%{
I = struct('A',0.3,...
    'B',0.01,...
    'm_dot',100,...
    'M_0',1.6,...
    'M_6',0.4,...
    'day','Standard',...
    'Alt',30000,...
    'Tt_4',3200,...
    'Tt_7',3600,...
    'h_PR',18000,...
    'CTOL',0.01,...
    'CTOH',0.00,...
    'e1',0.05,...
    'e2',0.05,...
    'e_f',0.890,... % Fan polytropic efficiency (e_f). [Dimensionless]
    'e_cL',0.890,... % LPC polytropic efficiency (e_cL). [Dimensionless]
    'e_cH',0.900,... % HPC polytropic efficiency (e_cH). [Dimensionless]
    'e_tH',0.890,... % HPT polytropic efficiency (e_tH). [Dimensionless]
    'e_tL',0.910,... % LPT polytropic efficiency (e_tL). [Dimensionless]
    'eta_b',0.98,... % Burner adiabatic efficiency (eta_b). [Dimensionless]
    'eta_AB',0.97,... % Afterburner adiabatic efficiency (eta_AB). [Dimensionless]
    'eta_mL',0.99,... % LP spool mechanical efficiency (eta_mL). [Dimensionless]
    'eta_mH',0.98,... % HP spool mechanical efficiency (eta_mH). [Dimensionless]
    'eta_mPH',0.98,... % HP PTO mechanical efficiency (eta_mPH). [Dimensionless]
    'eta_mPL',0.98,... % LP PTO mechanical efficiency (eta_mPL). [Dimensionless]
    'P0_P9',1,...
    'pi_AB',0.96,...
    'pi_f',3.5,...
    'pi_c',16,...
    'pi_cL',3.5,...
    'pi_b',0.97,...
    'pi_d_max',0.97,... % Maximum diffuser pressure ratio (pi_d_max). [Dimensionless]
    'pi_M_max',0.97,... % Maximum mixer pressure ratio (pi_M_max). [Dimensionless]
    'pi_n',0.98);
tic
R = ONX_LBTF_VSH(I)
toc
%}

M_0 = 1.1;
Alt = 30000;

tic
AB = 1.0;
O = OFX_LBTF_VH(R,I,M_0,Alt,AB)
toc
%}
function O = OFX_LBTF_VH(R,I,M_0,Alt,AB) %617
%Data Preallocations and Initial value setting.
O = struct('Alt',Alt,'A',0,'f',R.f,'ht_0',0,'ht_2',0,'ht_3',0,'ht_4',0,'ht_45',0,...
    'M_0',M_0,'M_6',R.M_6,'M_6A',R.M_6A,'Cp_6A',0,'gamma_6A',0,...
    'm_dot',R.m_dot,...
    'P_0',0,'T_0',0,'Tt_4',0,'Tt_45',0,'Tt_5',0,'Pt16_Pt6',0,'pi_ABdry',0,'pi_c',R.pi_c*1.1,'pi_cH',R.pi_cH,...
    'pi_cL',R.pi_cL,'pi_d',0,'pi_f',R.pi_f,'pi_tH',R.pi_tH,...
    'pi_tL',R.pi_tL,'pi_r',0,'tau_cH',R.tau_cH,'tau_cL',R.tau_cL,'tau_f',R.tau_f,...
    'tau_m1',R.tau_m1,'tau_m2',R.tau_m2,'tau_r',0,'tau_tH',R.tau_tH,...
    'tau_tL',R.tau_tL,'theta_0',0);

%Preliminary Computations
g_c = 32.174;
[O.T_0,~,O.P_0,rho_0] = ATMOSphere(O.Alt,'BE',I.day);
[~,h_0,Pr_0,Phi_0,Cp_0,R_0,gamma_0,a_0] = FAIR(1,O.T_0,0,'BE');
V_0 = O.M_0*a_0;
O.ht_0 = h_0 + (V_0^2)/(2*g_c*778.16);
O.tau_r = O.ht_0/h_0; %OFX Ram enthalpy ratio.
O.theta_0 = O.tau_r*O.T_0/518.67; %OFX Throttle Ratio.
[Tt_0,~,Prt_0,Phit_0,~,~,~,~] = FAIR(2,O.ht_0,0,'BE');
O.pi_r = Prt_0/Pr_0; %OFX Ram pressure ratio.
eta_R_spec = MIL_E_5008B(O.M_0);
O.pi_d = I.pi_d_max*eta_R_spec;
O.ht_2 = O.ht_0;
Prt_2 = Prt_0;
O.A = R.A*O.T_0*O.tau_r/(R.tau_r*R.T_0); %Initial guess for bypass.
O.pi_ABdry = 1-(1-I.pi_AB)/2;
if O.theta_0 >= R.theta_0
    O.Tt_4 = I.Tt_4;
    Tt_4_guess = I.Tt_4;
else
    Tt_4_guess = I.Tt_4*O.theta_0/R.theta_0;
    O.Tt_4 = Tt_4_guess;
end
[~,O.ht_4,Prt_4,Phit_4,~,~,~,~] = FAIR(1,O.Tt_4,O.f,'BE');
O.ht_45 = O.ht_4*O.tau_m1*O.tau_tH*O.tau_m2;
f_45 = O.f/(1+O.f+(R.e1+R.e2)/R.lmbme1me2);
[O.Tt_45,~,Prt_45,Phit_45,~,~,~,~] = FAIR(2,O.ht_45,f_45,'BE');
ht_5 = O.ht_45*O.tau_tL;
[O.Tt_5,~,Prt_5,Phit_5,~,~,~,~] = FAIR(2,ht_5,f_45,'BE');

%STAGE 1
k_max = 100;
k = 0;
k_throttle = 100;
while O.pi_c > R.pi_c && k_throttle > 0
    O.Tt_4 = Tt_4_guess*k_throttle/100;
    m_error = 1;
    while m_error > 0.001 && k <= k_max
        M6_error = 1;
        while M6_error > 0.0005 && k <= k_max
            a_error = 1;
            while a_error > 0.001 && k <= k_max
                req3 = true;
                while req3 == true && k <= k_max
                    O.ht_3 = O.ht_0*O.tau_cL*O.tau_cH;
                    [Tt_3,~,Prt_3,Phit_3,~,~,~,~] = FAIR(2,O.ht_3,0,'BE');
                    O.A_prime = O.A/((1+O.f)*R.lmbme1me2+R.e1+R.e2);%%%%%%%%%%
                    [~,O.ht_4,Prt_4,Phit_4,~,~,~,~] = FAIR(1,O.Tt_4,O.f,'BE');
                    [O.pi_tH,O.tau_tH,O.Tt_45] = TURBC(O.Tt_4,O.ht_4,O.f,R.A4_A45,R.M_4,R.M_45,R.eta_tH,O.Tt_45,O.ht_3,R.lmbme1me2,R.e1,R.e2);
                    [O.pi_tL,O.tau_tL,O.Tt_5] = TURB(O.Tt_45,f_45,R.A45_A6,R.M_45,O.M_6,R.eta_tL,O.Tt_5);
                    [~,ht_5,Prt_5,Phit_5,~,~,~,~] = FAIR(1,O.Tt_5,f_45,'BE');
                    O.tau_lambda = O.ht_4/h_0;
                    O.tau_f = 1 + ((1-O.tau_tL)*I.eta_mL*(R.lmbme1me2*(1+O.f)*O.tau_lambda*O.tau_tH/O.tau_r+(R.e1*O.tau_tH+R.e2)*O.tau_cL*O.tau_cH)-(1+O.A)*R.P_TOL/(O.tau_r*I.eta_mPL*O.m_dot*h_0))/((R.tau_cL-1)/(R.tau_f-1)+O.A);
                    O.tau_cL = 1 + (O.tau_f-1)*(R.tau_cL-1)/(R.tau_f-1);
                    O.tau_cH = (1 + (1-O.tau_tH)*I.eta_mH*(R.lmbme1me2*(1+O.f)*O.tau_lambda/(O.tau_r*O.tau_cL)) - (1+O.A)*R.P_TOH/(O.tau_r*O.tau_cL*I.eta_mPH*O.m_dot*h_0))/(1-R.e1*(1-O.tau_tH)*I.eta_mH);
                    ht_13 = O.ht_2*O.tau_f;
                    ht_25 = O.ht_2*O.tau_cL;
                    O.ht_3 = ht_25*O.tau_cH;
                    ht_13i = O.ht_2*(1+R.eta_f*(O.tau_f-1));
                    ht_25i = O.ht_2*(1+R.eta_cL*(O.tau_cL-1));
                    ht_3i = ht_25*(1+R.eta_cH*(O.tau_cH-1));
                    [Tt_13,~,Prt_13,Phit_13,~,~,~,~] = FAIR(2,ht_13,0,'BE');
                    [Tt_25,~,Prt_25,Phit_25,~,~,~,~] = FAIR(2,ht_25,0,'BE');
                    [Tt_3,~,Prt_3,Phit_3,~,~,~,~] = FAIR(2,O.ht_3,0,'BE');
                    [Tt_13i,~,Prt_13i,Phit_13i,~,~,~,~] = FAIR(2,ht_13i,0,'BE');
                    [Tt_25i,~,Prt_25i,Phit_25i,~,~,~,~] = FAIR(2,ht_25i,0,'BE');
                    [Tt_3i,~,Prt_3i,Phit_3i,~,~,~,~] = FAIR(2,ht_3i,0,'BE');
                    O.pi_f = Prt_13i/Prt_2;
                    O.pi_cL = Prt_25i/Prt_2;
                    O.pi_cH = Prt_3i/Prt_25;
                    O.pi_c = O.pi_cL*O.pi_cH;
                    O.f = FUEL_NR(O.ht_3,I.eta_b,I.h_PR,O.Tt_4);
                    O.tau_m1 = (R.lmbme1me2*(1+O.f)+R.e1*O.tau_r*O.tau_cL*O.tau_cH/O.tau_lambda)/(R.lmbme1me2*(1+O.f)+R.e1);
                    O.tau_m2 = (R.lmbme1me2*(1+O.f)+R.e1+R.e2*(O.tau_r*O.tau_cL*O.tau_cH/(O.tau_lambda*O.tau_m1*O.tau_tH)))/(R.lmbme1me2*(1+O.f)+R.e1+R.e2);
                    ht_6 = ht_5;
                    Tt_6 = O.Tt_5;
                    ht_16 = ht_13;
                    Tt_16 = Tt_13;
                    Pt_6 = O.P_0*O.pi_r*O.pi_d*O.pi_cL*O.pi_cH*I.pi_b*O.pi_tH*O.pi_tL;
                    f_45 = O.f/(1+O.f+(R.e1+R.e2)/R.lmbme1me2);
                    [TtT6,PtP6,O.MFP_6] = MASSFP_NR(Tt_6,f_45,O.M_6);
                    P_6 = Pt_6/PtP6;
                    T_6 = Tt_6/TtT6;
                    Pt_16 = O.P_0*O.pi_r*O.pi_d*O.pi_f;
                    PtP16 = Pt_16/P_6; %Kutta condition
                    [~,PtP,~] = MASSFP_NR(Tt_16,0,1); %PtP at M = 1
                    if PtP16 > PtP || PtP16 < 1
                        O.M_6 = O.M_6-0.01;
                        k = k+1;
                    else
                        req3 = false;
                    end
                    %fprintf('loop1\n')
                end
                [O.M_16,TtT16,~,MFP_16] = RGCOMPR_NR(3,Tt_16,0,PtP16);
                T_16 = Tt_16/TtT16;
                A_prime = R.A16_A6*Pt_16*MFP_16*sqrt(Tt_6/Tt_16)/(Pt_6*O.MFP_6);
                O.Pt16_Pt6 = Pt_16/Pt_6;
                a_error = abs((A_prime-O.A_prime)/O.A_prime);
                O.A = A_prime*((1+O.f)*R.lmbme1me2+R.e1+R.e2);
                k = k+1;
                %fprintf('loop2\n')
            end
            [~,h_6,Pr_6,Phi_6,Cp_6,R_6,gamma_6,a_6] = FAIR(1,T_6,f_45,'BE');
            [~,h_16,Pr_16,Phi_16,Cp_16,R_16,gamma_16,a_16] = FAIR(1,T_16,0,'BE');
            ht_6A = (ht_6+O.A_prime*ht_16)/(1+O.A_prime);
            O.tau_M = ht_6A/ht_6;
            O.f_6A = f_45/(1+O.A_prime);
            [O.M_6A,MFP_6A,Tt_6A,O.gamma_6A,O.Cp_6A] = MIXER_VSH(O.M_6,O.M_16,gamma_16,gamma_6,R_6,O.A_prime,R.A16_A6,T_6,O.f_6A,ht_6A);
            O.pi_M_ideal = sqrt(Tt_6A/Tt_6)*O.MFP_6*(1+O.A_prime)/(MFP_6A*(1+R.A16_A6));
            O.pi_M = O.pi_M_ideal*I.pi_M_max;
            PtP9 = O.pi_r*O.pi_d*O.pi_cL*O.pi_cH*I.pi_b*O.pi_tH*O.pi_tL*O.pi_M*O.pi_ABdry*I.pi_n*I.P0_P9;
            [O.M_9,TtT9,~,MFP_9] = RGCOMPR_NR(3,Tt_6A,O.f_6A,PtP9);
            if O.M_9 > 1
                O.M_8 = 1;
            else
                O.M_8 = O.M_9;
            end
            [TtT_8,PtP_8,MFP_8] = MASSFP_NR(Tt_6A,O.f_6A,O.M_8);
            O.MFP_6 =  MFP_8*O.pi_M*O.pi_ABdry*R.A8_A6_dry*sqrt(Tt_6/Tt_6A)/(1+O.A_prime);
            [M_6new,TtT_6,PtP_6,~] = RGCOMPR_NR(5,Tt_6,f_45,O.MFP_6);
            M6_error = abs(M_6new-O.M_6);
            if M6_error > 0.0005
                if O.M_6 > M_6new
                    O.M_6 = O.M_6 - 0.0001;
                else
                    O.M_6 = O.M_6 + 0.002;
                end
            end
            %fprintf('loop3\n')
        end
        [~,TtT,PtP,MFP_4] = RGCOMPR_NR(1,O.Tt_4,O.f,1);
        m_new = R.m_dot*(1+R.f)*O.P_0*(1+O.A)*O.pi_r*O.pi_d*O.pi_c*MFP_4*sqrt(R.Tt_4/O.Tt_4)/((1+O.f)*(R.P_0*(1+R.A)*R.pi_r*R.pi_d*R.pi_c*R.MFP_4));
        m_error = abs((m_new-O.m_dot)/O.m_dot);
        O.m_dot = m_new;
        k = k + 1;
        %fprintf('loop4\n')
    end
    O.Tt_7 = AB*(R.Tt_7-Tt_6A) + Tt_6A;
    if AB == 0
        O.f_AB =0;
    else
        O.f_AB = FUEL_NR(ht_6A,I.eta_b,I.h_PR,O.Tt_7);
    end
    O.f_7 = O.f_6A + O.f_AB;
    O.pi_AB = O.pi_ABdry + 0.01*AB*(I.pi_AB-O.pi_ABdry);
    PtP9 = O.pi_r*O.pi_d*O.pi_cL*O.pi_cH*I.pi_b*O.pi_tH*O.pi_tL*O.pi_M*O.pi_AB*I.pi_n*I.P0_P9;
    Tt_9 = O.Tt_7;
    [M_9,TtT9,~,MFP_9] = RGCOMPR_NR(3,Tt_9,O.f_7,PtP9);
    m_dot_9 = O.m_dot*(1+O.f_7)*(1-R.B/(1+O.A));
    Pt_9 = O.P_0*O.pi_r*O.pi_d*O.pi_cL*O.pi_cH*I.pi_b*O.pi_tH*O.pi_tL*O.pi_M*O.pi_AB*I.pi_n;
    O.A_9 = m_dot_9*sqrt(Tt_9)/(Pt_9*MFP_9);
    T_9 = Tt_9/TtT9;
    [~,h_9,Pr_9,Phi_9,Cp_9,R_9,gamma_9,a_9] = FAIR(1,T_9,O.f_7,'BE');
    V_9 = M_9*a_9;
    O.f_0 = (O.f*R.lmbme1me2 + O.f_AB*(1+O.A-R.B))/(1+O.A);
    O.F_m_dot = (a_0/g_c)*((1+O.f_0 - I.B/(1+I.A))*(V_9/a_0 + R_9*T_9*a_0*(1-I.P0_P9)/(R_0*R.T_0*V_9*gamma_0)) - I.M_0); % Uninstalled Specific Thrust (AB = on). [lbf/lbm*s]
    O.S = 3600*O.f_0/O.F_m_dot;
    O.F = O.m_dot*O.F_m_dot;
    [~,TtT0,PtP0,MFP_0] = RGCOMPR_NR(1,Tt_0,0,O.M_0);
    O.A_0 = O.m_dot*sqrt(Tt_0)/(O.P_0*PtP0*MFP_0);
    O.LP_RPM = 100*sqrt(h_0*O.tau_r*(O.tau_f-1)/(R.h_0*R.tau_r*(R.tau_f-1)));
    O.HP_RPM = 100*sqrt(h_0*O.tau_r*O.tau_cL*(O.tau_cH-1)/(R.h_0*R.tau_r*R.tau_cL*(R.tau_cH-1)));
    O.eta_P = (2*g_c*O.M_0*O.F_m_dot/a_0)/((1 + O.f_0 - R.B/(1+O.A))*((V_9/a_0)^2)-O.M_0^2);
    O.eta_TH =(((1+O.f_0-I.B/(1+I.A))*(V_9^2) - V_0^2)/(2*g_c*778.16)+(I.CTOL+I.CTOH)*h_0)/(O.f_0*I.h_PR);
    k_throttle = k_throttle - 1;
end
    fprintf('success!\n')
    fprintf('%5f',k_throttle)
end
%}