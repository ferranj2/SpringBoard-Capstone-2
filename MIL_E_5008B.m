function eta_R_spec = MIL_E_5008B(M_0)
if M_0 <= 1
    eta_R_spec = 1;
elseif 1 < M_0 && M_0 <5
    eta_R_spec = 1 - 0.075*power(M_0 - 1,1.35);
else
    eta_R_spec = 800/(M_0^4 +935);
end
end