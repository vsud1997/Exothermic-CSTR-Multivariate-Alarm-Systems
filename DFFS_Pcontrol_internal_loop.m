function [Cs_1_saved, Ts_1_saved, Tc_1_saved, Fc_1_saved, N, probability_N_i] = DFFS_Pcontrol_internal_loop(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Fc_ub, lambda_i, N_i, N_traj, traj_length, Ts_saved,Cs_saved,Fc_saved,Tc_saved)

A = 1/tau; % 1/ResidenceTime =  1/0.5 = 2
ko = 17.038; % Rate constant
E = 1.50e+04; % Actication energy
R = 8.314; % Gas Constant
alpha = 0.075; 
beta = 9370.9;
h = 0.01;

Cs_1_saved = [];
Ts_1_saved = [];
dW_1_saved = [];
ts_1_saved = [];
Fc_1_saved = [];
Tc_1_saved = [];

N = 0;
b = 1/tau;
b1 = 1;

pd = makedist('Normal',0,sqrt(var)); 

for k = 1:N_traj % Trajetories for crossing lambda_1 K

   rc = randi(N_i); %Random crossing selected
   Ts_init = Ts_saved(rc); % The rc^(th) element will be the initial values
   Cs_init = Cs_saved(rc);
   Tc_init = Tc_saved(rc);
   Fc_init = Fc_saved(rc);
   ts_init = 0;

   Cs = [];
   Ts = [];
   ts = [];
   Fc = [];
   Tc = [];
  
  Cs(1) = Cs_init;
  Ts(1) = Ts_init;
  Fc(1) = Fc_init;
  Tc(1) = Tc_init;
  ts(1) = ts_init;

    C = Cs_init;
    T = Ts_init;
    Tc = Tc_init;
    Fc = Fc_init;
    eP_ = b1*Tsp - T;
    F = Fc_init;

    for i = 1:traj_length-1
    % Initiating trajectories from each randomly selected configuration
    % points. Assumed 500 steps.
    C = Cs(i) ;
    T = Ts(i) ;   
    Tc_s = Tc(i);
    eP = b*Tsp - T;
    %F = (Fc_0 + Kc*(T - Tsp));
    F = F - (Kc*(eP - eP_));
    dW = random(pd);
    noise = dW ;
  
    %% CW constraints should come here!
    F = max(Fc_lb, min(F, Fc_ub));

    %% Create functions at current time-step
    a = A*(2-C) - ko*C*exp(-E/R/T);
    a_T = A*(300-T) + alpha*(Tc_s-T) + beta*C*exp(-E/R/T); % RHS of energy balance equation to iterate
    a_Tc = F*(300-Tc_s)/10 - alpha*(Tc_s - T);
    %% Update for next time-step
    
    Cs(i+1) = C + a*h + b*dW*sqrt(h); % Updating concentration for next time step
    Ts(i+1) = T + a_T*h  ; % Updating temperature for next time step
    Tc(i+1) = Tc_s + a_Tc*h;
    ts(i+1) = ts(i) + h;
    Fc(i+1) = F;

    eP_ = eP;



      if  Ts(i+1)< Ts(i) && Ts(i+1)<lambda_i && Ts(i)>lambda_i 
         N = N + 1;
         Cs_1_saved = [Cs_1_saved, Cs(i+1)];
         Ts_1_saved = [Ts_1_saved, Ts(i+1)];
         dW_1_saved = [dW_1_saved, noise];
         ts_1_saved = [ts_1_saved, ts(i+1)];
         Fc_1_saved = [Fc_1_saved, Fc(i+1)];
         Tc_1_saved = [Tc_1_saved, Tc(i+1)];
         break
      end
    end
end
probability_N_i = N/N_traj;

end