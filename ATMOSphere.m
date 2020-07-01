%TO-DO:
%1) Code-in the cold day model.
%2) Code-in the tropical day model.

%----------------------------------------
%AE 435: Air-Breathing Preliminary Design
%RFP: AEDsys "ATMOS" Code
%By Jesus Ferrand
%Submitted to: Mark Ricklick Phd.
%----------------------------------------

%NOTES:
%AEDsys' ATMOS program only has data up to h = 30km (98.4252 kft) for its
%standard day model. The data for the cold, hot, and tropic days is limited
%to heights of 22.5km (73.8189 kft), 20.5km (67.2572 kft), and 21km
%(68.8976 kft) respectively. This code will accept inputs outside these
%height limits but bear in mind that the output is extrapolated.

%CODE:
%==========================================================================
%[T,a,P,rho] = ATMOSphere(10000,'SI','Standard');
%[T,a,P,rho] = ATMOSphere(10000/0.3048,'BE','Standard');
%[T,a,P,rho] = ATMOSphere(10000,'SI','Standard','on')

function [T,a,P,rho] = ATMOSphere(h,units,day)
if strcmpi(units,'BE') == 1
    h = h*0.3048; % Input height in BE units converted to SI units. [m]
end
r_0 = 6356.577*1000; % Earth's radius. [m]
g_0 = 9.80665; % Earth's acceleration due to gravity at its surface. [m/s^2]
R = 8.31432; % Universal gas constant. [J/mol-K]
W_0 = 28.9644/1000; % Air molecular weight. [kg/mol]
P_0 = 101325; % SLS pressure. [Pa]
z = h*r_0/(r_0+h); % "Height Parameter." [m]

if strcmpi(day,'Standard') == 1
    T_0 = 288.15; % Surface Temperature. [K]
    c = [0*1000,-6.5/1000;...
        11*1000,0/1000;...
        20*1000,1/1000;...
        32*1000,2.8/1000;...
        47*1000,0/1000;...
        51*1000,-2.8/1000;...
        71*1000,-2/1000;...
        84.582*1000,0/1000];
    
    T_1 = T_0 + c(1,2)*(c(2,1)-c(1,1)); % [K]
    T_2 = T_1 + c(2,2)*(c(3,1)-c(2,1)); % [K]
    T_3 = T_2 + c(3,2)*(c(4,1)-c(3,1)); % [K]
    T_4 = T_3 + c(4,2)*(c(5,1)-c(4,1)); % [K]
    T_5 = T_4 + c(5,2)*(c(6,1)-c(5,1)); % [K]
    T_6 = T_5 + c(6,2)*(c(7,1)-c(6,1)); % [K]
    
    P_1 = P_0*power(T_0/T_1,g_0*W_0/(R*c(1,2))); % [Pa]
    P_2 = P_1*exp(-g_0*W_0*(c(3,1)-c(2,1))/(R*T_1)); % [Pa]
    P_3 = P_2*power(T_2/T_3,g_0*W_0/(R*c(3,2))); % [Pa]
    P_4 = P_3*power(T_3/T_4,g_0*W_0/(R*c(4,2))); % [Pa]
    P_5 = P_4*exp(-g_0*W_0*(c(6,1)-c(5,1))/(R*T_4)); % [Pa]
    P_6 = P_5*power(T_5/T_6,g_0*W_0/(R*c(6,2))); % [Pa]
    %P_7 = P_6*power(T_6/T_7,g_0*W_0/(R*c(7,2))); % [Pa]
    
    if z >= 0 && z < c(2,1)
        T_i = T_0;
        P_i = P_0;
        L_i = c(1,2);
        z_i = c(1,1);
    elseif z >= c(2,1) && z < c(3,1)
        T_i = T_1;
        P_i = P_1;
        L_i = c(2,2);
        z_i = c(2,1);
    elseif z >= c(3,1) && z < c(4,1)
        T_i = T_2;
        P_i = P_2;
        L_i = c(3,2);
        z_i = c(3,1);
    elseif z >= c(4,1) && z < c(5,1)
        T_i = T_3;
        P_i = P_3;
        L_i = c(4,2);
        z_i = c(4,1);
    elseif z >= c(5,1) && z < c(6,1)
        T_i = T_4;
        P_i = P_4;
        L_i = c(5,2);
        z_i = c(5,1);
    elseif z >= c(6,1) && z < c(7,1)
        T_i = T_5;
        P_i = P_5;
        L_i = c(6,2);
        z_i = c(6,1);
    elseif z >= c(7,1) && z < c(8,1)
        T_i = T_6;
        P_i = P_6;
        L_i = c(7,2);
        z_i = c(7,1);
    end
    T = T_i + L_i*(z-z_i); % Temperature @ input height. [K]
    P = pressure_lapse(L_i,T,T_i,P_i,z,z_i,R,g_0,W_0); % Static Pressure @ input height. [Pa]
    rho = P*W_0/(R*T); % Static Density @ input height. [kg/m^3]
    a = sqrt(1.4*R*T/W_0); % Speed of Sound @ input height. [m/s]
    
elseif strcmpi(day,'Hot') == 1
    T_0 = 312.60; % Surface Temperature. [K]
    c = [0*1000,-7/1000;...
        12*1000,0.8/1000;...
        20.5*1000,1.4/1000];
    T_1 = T_0 + c(1,2)*(c(2,1)-c(1,1)); % [K]
    T_2 = T_1 + c(2,2)*(c(3,1)-c(2,1)); % [K]
    P_1 = P_0*power(T_0/T_1,g_0*W_0/(R*c(1,2))); % [Pa]
    P_2 = P_1*power(T_1/T_2,g_0*W_0/(R*c(2,2))); % [Pa]
    if h >= 0 && h < c(2,1) %0km <= h < 12km
        T_i = T_0;
        P_i = P_0;
        L_i = c(1,2);
        h_i = c(1,1);
    elseif h >= c(2,1) && h < c(3,1) %12km <= h < 20.5km
        T_i = T_1;
        P_i = P_1;
        L_i = c(2,2);
        h_i = c(2,1);
    elseif h >= c(3,1) %20.5km <= h
        T_i = T_2;
        P_i = P_2;
        L_i = c(3,2);
        h_i = c(3,1);
    end
    z_i = h_i*r_0/(r_0+h_i);
    T = T_i + L_i*(h-h_i); % Temperature @ input height. [K]
    P = pressure_lapse(L_i,T,T_i,P_i,z,z_i,R,g_0,W_0); % Static Pressure @ input height. [Pa]
    rho = P*W_0/(R*T); % Static Density @ input height. [kg/m^3]
    a = sqrt(1.4*R*T/W_0); % Speed of Sound @ input height. [m/s]
elseif strcmpi(day,'Cold') == 1
    T_0 = 222.1; % Surface Temperature. [K]
elseif strcmpi(day,'Tropical') == 1
    T_0 = 305.27; % Surface Temperature. [K]
end
if strcmpi(units,'BE') == 1
    T = T*1.8; % Static Temperature @ input height. [R]
    a = a/0.3048; % Speed of sound @ input height. [ft/s]
    P = P*2116.2/101325; % Static Pressure @ input height. [psf]
    rho = rho/((100^3)/((12*2.54)^3))/(0.45359*32.174); % Static Density @ input height. [slug/ft^3]
end
    function P = pressure_lapse(L_i,T,T_i,P_i,z,z_i,R,g_0,W_0)
        if L_i == 0
            P = P_i*exp(-g_0*W_0*(z-z_i)/(R*T_i));
        else
            P = P_i*power(T_i/T,g_0*W_0/(R*L_i));
        end
    end
end