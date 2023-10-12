function [A,B] = jacobian_CSTR_P_controlled(k0, tau, Kc, Fc_0, Tc_0,R, alpha, beta, V_j, E, tBounds, C_ss, T_ss, Fc_ss, Tc_ss)

A = zeros(3,3);
B = zeros(3,1);

A(1,1) = -1/tau - k0*exp(-E/R/T_ss);
A(1,2) = -k0*C_ss*(E/R/T_ss^2)*exp(-E/R/T_ss);
A(1,3) = 0;

A(2,1) = beta*exp(-E/R/T_ss);
A(2,2) = -1/tau - alpha + beta*C_ss*(E/R/T_ss^2)*exp(-E/R/T_ss);
A(2,3) = alpha;

A(3,1) = 0;
A(3,2) = alpha;
A(3,3) = - alpha - (Fc_ss/V_j);

B(1) = 0;
B(2) = 0;
B(3) = (Tc_0 - Tc_ss)/V_j;

end