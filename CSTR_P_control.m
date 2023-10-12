function [Cs, Ts, Fc, Tc] = CSTR_P_control(tau, tBounds, init_var, Fc_0, Tc_0, Kc)


A = 1/tau; % 1/ResidenceTime =  1/0.5 = 2
ko = 17.038; % Rate constant
E = 1.50e+04; % Actication energy
R = 8.314; % Gas Constant
alpha = 0.075; 
beta = 9370.9;
h = 0.01;
Tsp = 800;
N_steps = (tBounds(2) - tBounds(1))/h + 1;
Cs  = zeros(N_steps,1);
Ts = zeros(N_steps,1);
Fc = zeros(N_steps,1);
Tc = zeros(N_steps,1);

Cs(1) = init_var(1);
Ts(1) = init_var(2);
Fc(1) = init_var(3);
Tc(1) = init_var(4);

for i = 2:N_steps

    C = Cs(i-1) ;
    T = Ts(i-1) ;
    a = A*(2-C) - ko*C*exp(-E/R/T);
    Tc_s = Tc(i-1);
    
    
    a_T = A*(300-T) + alpha*(Tc_s-T) + beta*C*exp(-E/R/T); % RHS of energy balance equation to iterate

    F = (Fc_0 + Kc*( T - Tsp));
    
    a_Tc = F*(300-Tc_s)/10 - alpha*(Tc_s - T);
    
    Cs(i) = C + a*h ; % Updating concentration for next time step
    
    Ts(i) = T + a_T*h  ; % Updating temperature for next time step
    
    Tc(i) = Tc_s + a_Tc*h;

  
    Fc(i) = F;
    

end

end