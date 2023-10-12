function [Ts_new_guess, Cs_new_guess, Fc_new_guess, Tc_new_guess, noise_new_guess,ts_new_guess,overall_prob,N1, Ts_1_saved,Cs_1_saved,dW_1_saved, Fc_1_saved, Tc_1_saved, ts_1_saved, cell_Ts_1,cell_Cs_1, cell_Fc_1, cell_Tc_1,  cell_noise_1, cell_N1,cell_ts_1,cell_Ts_2,cell_Cs_2,cell_Fc_2, cell_Tc_2, cell_noise_2,cell_N2,cell_ts_2,cell_Ts_3,cell_Cs_3, cell_Fc_3, cell_Tc_3, cell_noise_3, cell_N3,cell_ts_3, cell_Ts_4,cell_Cs_4, cell_Fc_4, cell_Tc_4, cell_noise_4, cell_N4,cell_ts_4, cell_Ts_5,cell_Cs_5, cell_Fc_5, cell_Tc_5, cell_noise_5, cell_N5,cell_ts_5, cell_Ts_6,cell_Cs_6, cell_Fc_6, cell_Tc_6, cell_noise_6, cell_N6,cell_ts_6, cell_Ts_7,cell_Cs_7, cell_Fc_7, cell_Tc_7, cell_noise_7, cell_N7,cell_ts_7] = outer_BG_FFS_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Cs_saved, Ts_saved, Fc_saved, Tc_saved, noise_saved, ts_saved, N_0, k0, k1, k2, k3 , k4 , k5 ,k6, k7)

A = 1/tau; % 1/ResidenceTime =  1/0.5 = 2
ko = 17.038; % Rate constant
E = 1.50e+04; % Actication energy
R = 8.314; % Gas Constant
alpha = 0.075; 
beta = 9370.9;
h = 0.01;
N_interface = 8;
transition_time = 16; %minutes
time_traj = 3*transition_time/N_interface;
traj_length = time_traj/h + 1;
b = 1/tau;
pd = makedist('Normal',0,sqrt(var)); 

lambda_1 = 620;
lambda_2 = 580;
lambda_3 = 548;
lambda_4 = 524;
lambda_5 = 492;
lambda_6 = 455;
lambda_7 = 425;
lambda_8 = 400;

N1 = 0; 
Ts_1_saved = [];
Cs_1_saved = [];
dW_1_saved = [];
ts_1_saved = [];
Fc_1_saved = [];
Tc_1_saved = [];


rc_1 = randi(N_0);

Ts_new_guess = Ts_saved(rc_1); 
Cs_new_guess = Cs_saved(rc_1);
noise_new_guess = noise_saved(rc_1);
%ts_new_guess = ts_saved(rc_1);
ts_new_guess = 0;
Fc_new_guess = Fc_saved(rc_1);
Tc_new_guess = Tc_saved(rc_1);

tic

for j = 1:k0 % Trajetories for crossing at lambda_0
    
   Cs  = [];
   Ts = [];
   ts = [];
   Fc = [];
   Tc = [];
  
   Cs(1) = Cs_new_guess;
   Ts(1) = Ts_new_guess;
   ts(1) = ts_new_guess;
   Fc(1) = Fc_new_guess;
   Tc(1) = Tc_new_guess;

    for i = 2:traj_length
    % Initiating trajectories from each randomly selected configuration
    % points. Assumed 500 steps.
    C = Cs(i-1) ;
    T = Ts(i-1) ;
    Tc_s = Tc(i-1);
    F = (Fc_0 + Kc*(T - Tsp));
    %F = Fc(i-1);
    dW = random(pd);
    noise = dW ;

    %% CW constraints should come here!
    if F < Fc_lb
        F = Fc_lb;
    end
    %% Update for next time-step
    a = A*(2-C) - ko*C*exp(-E/R/T);
    a_T = A*(300-T) + alpha*(Tc_s-T) + beta*C*exp(-E/R/T); % RHS of energy balance equation to iterate
    a_Tc = F*(300-Tc_s)/10 - alpha*(Tc_s - T);
    
    %% Update for next time-step
    Cs(i) = C + a*h + b*dW*sqrt(h); % Updating concentration for next time step

    Ts(i) = T + a_T*h  ; % Updating temperature for next time step
    
    Tc(i) = Tc_s + a_Tc*h;

    ts(i) = ts(i-1) + h;

    Fc(i) = F;
    
    if  Ts(i)< Ts(i-1) && Ts(i)<lambda_1 && Ts(i-1)>lambda_1 %&& Cs_1(i) > solving_for_conc_as_func_of_T(tau, 580)- 2*sqrt(var) && Cs_1(i) < solving_for_conc_as_func_of_T(tau, 580)+ 2*sqrt(var)% If The trajectory crosses over the first order parameter of 580 K
         
         N1 = N1 + 1;
         Cs_1_saved = [Cs_1_saved, Cs(i)];
         Ts_1_saved = [Ts_1_saved, Ts(i)];
         dW_1_saved = [dW_1_saved, noise];
         ts_1_saved = [ts_1_saved, ts(i)];
         Fc_1_saved = [Fc_1_saved, Fc(i)];
         Tc_1_saved = [Tc_1_saved, Tc(i)];

         break
  
     end
     
    end

   %ts_prior_cell_k0{j} = ts_prior;
   %Ts_prior_cell_k0{j} = Ts_prior;
   %Cs_prior_cell_k0{j} = Cs_prior;
   %Fc_prior_cell_k0{j} = Fc_prior;
   %Tc_prior_cell_k0{j} = Tc_prior;
    
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

[cell_Ts_1, cell_Cs_1,cell_noise_1,cell_N1,cell_ts_1,cell_Fc_1,cell_Tc_1, Temp_1_saved,Conc_1_saved,sum_N1,noise_1_saved,time_1_saved,Flow_1_saved,Temp_CW_1_saved] = inner_BGFFS_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0,Fc_lb,lambda_2,N1,k1,Ts_1_saved,Cs_1_saved,Fc_1_saved,Tc_1_saved,ts_1_saved);

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

[cell_Ts_2, cell_Cs_2,cell_noise_2,cell_N2,cell_ts_2,cell_Fc_2,cell_Tc_2, Temp_2_saved,Conc_2_saved,sum_N2,noise_2_saved,time_2_saved,Flow_2_saved,Temp_CW_2_saved] = inner_BGFFS_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0,Fc_lb,lambda_3,sum_N1,k2,Temp_1_saved,Conc_1_saved,Flow_1_saved,Temp_CW_1_saved, time_1_saved);

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

[cell_Ts_3, cell_Cs_3,cell_noise_3,cell_N3,cell_ts_3,cell_Fc_3,cell_Tc_3, Temp_3_saved,Conc_3_saved,sum_N3,noise_3_saved,time_3_saved,Flow_3_saved,Temp_CW_3_saved] = inner_BGFFS_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0,Fc_lb,lambda_4,sum_N2,k3,Temp_2_saved,Conc_2_saved,Flow_2_saved,Temp_CW_2_saved,time_2_saved);

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

[cell_Ts_4, cell_Cs_4,cell_noise_4,cell_N4,cell_ts_4,cell_Fc_4,cell_Tc_4, Temp_4_saved,Conc_4_saved,sum_N4,noise_4_saved,time_4_saved,Flow_4_saved,Temp_CW_4_saved] = inner_BGFFS_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0,Fc_lb,lambda_5,sum_N3,k4,Temp_3_saved,Conc_3_saved,Flow_3_saved,Temp_CW_3_saved,time_3_saved);

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

[cell_Ts_5, cell_Cs_5,cell_noise_5,cell_N5,cell_ts_5,cell_Fc_5,cell_Tc_5, Temp_5_saved,Conc_5_saved,sum_N5,noise_5_saved,time_5_saved,Flow_5_saved,Temp_CW_5_saved] = inner_BGFFS_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0,Fc_lb,lambda_6,sum_N4,k5,Temp_4_saved,Conc_4_saved,Flow_4_saved,Temp_CW_4_saved,time_4_saved);

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

[cell_Ts_6, cell_Cs_6,cell_noise_6,cell_N6,cell_ts_6,cell_Fc_6,cell_Tc_6, Temp_6_saved,Conc_6_saved,sum_N6,noise_6_saved,time_6_saved,Flow_6_saved,Temp_CW_6_saved] = inner_BGFFS_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0,Fc_lb,lambda_7,sum_N5,k6,Temp_5_saved,Conc_5_saved,Flow_5_saved,Temp_CW_5_saved,time_5_saved);

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

[cell_Ts_7, cell_Cs_7,cell_noise_7,cell_N7,cell_ts_7,cell_Fc_7,cell_Tc_7, Temp_7_saved,Conc_7_saved,sum_N7,noise_7_saved,time_7_saved,Flow_7_saved,Temp_CW_7_saved] = inner_BGFFS_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0,Fc_lb,lambda_8,sum_N6,k7,Temp_6_saved,Conc_6_saved,Flow_6_saved,Temp_CW_6_saved,time_6_saved);

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%overall_prob = sum_N6/(k0*k1*k2*k3*k4*k5*k6)

overall_prob = sum_N7/(k0*k1*k2*k3*k4*k5*k6*k7)

end

%[Ts_new_guess, Cs_new_guess, Fc_new_guess, Tc_new_guess, noise_new_guess,ts_new_guess,overall_prob,N1, Ts_1_saved,Cs_1_saved,dW_1_saved, Fc_1_saved, Tc_1_saved, ts_1_saved, cell_Ts_1,cell_Cs_1, cell_Fc_1, cell_Tc_1,  cell_noise_1, cell_N1,cell_ts_1,cell_Ts_2,cell_Cs_2,cell_Fc_2, cell_Tc_2, cell_noise_2,cell_N2,cell_ts_2,cell_Ts_3,cell_Cs_3, cell_Fc_3, cell_Tc_3, cell_noise_3, cell_N3,cell_ts_3, cell_Ts_4,cell_Cs_4, cell_Fc_4, cell_Tc_4, cell_noise_4, cell_N4,cell_ts_4, cell_Ts_5,cell_Cs_5, cell_Fc_5, cell_Tc_5, cell_noise_5, cell_N5,cell_ts_5, cell_Ts_6,cell_Cs_6, cell_Fc_6, cell_Tc_6, cell_noise_6, cell_N6,cell_ts_6, cell_Ts_7,cell_Cs_7, cell_Fc_7, cell_Tc_7, cell_noise_7, cell_N7,cell_ts_7] = outer_BG_FFS_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0, Cs_saved, Ts_saved, Fc_saved, Tc_saved, noise_saved, ts_saved, N_0, k0, k1, k2, k3 , k4 , k5 ,k6, k7)
