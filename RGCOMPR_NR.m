function [M,TtT,PtP,MFP] = RGCOMPR_NR(mode,Tt,f,IN)
switch mode
    case 1 %Mach Number known
        M = IN;
        [TtT,PtP,MFP] = MASSFP_NR(Tt,f,M);
    case {2,3}
        [~,ht,Prt,Phit,~,~,~,~] = FAIR(1,Tt,f,'BE');
        if mode == 2 %TtT known;
            TtT = IN;
            T = Tt/TtT;
            [~,h,Pr,Phi,Cp,R,gamma,a] = FAIR(1,T,f,'BE');
            PtP = Prt/Pr;
        else %PtP known;
            PtP = IN;
            Pr = Prt/PtP;
            [T,h,~,Phi,Cp,R,gamma,a] = FAIR(3,Pr,f,'BE');
            TtT = Tt/T;
        end
        V2 = 2*31.174*778.16*(ht-h);
        if V2 < 0
            M = 0;
            TtT = 1;
        else
            M = sqrt(V2)/a;
        end
        MFP = (M/PtP)*sqrt(gamma*TtT*32.174/R);
    case {4,5} %MFP Known
        MFP = IN;
        [TtT,PtP,M] = INV_MASSP_NR(mode,Tt,f,MFP);
end
end