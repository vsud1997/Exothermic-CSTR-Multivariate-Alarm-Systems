function [cell_Ts_1, cell_Cs_1,cell_noise_1,cell_N1,cell_ts_1,cell_Fc_1,cell_Tc_1, Temp_1_saved,Conc_1_saved,sum_N1,noise_1_saved,time_1_saved,Flow_1_saved,Temp_CW_1_saved] = inner_BGFFS_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0,Fc_lb, lambda_i,N_i,k_i,Ts_saved,Cs_saved,Fc_saved,Tc_saved,ts_saved)

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

cell_Ts_1 = cell(1, N_i);
cell_Cs_1 = cell(1, N_i);
cell_noise_1 = cell(1, N_i);
cell_N1 = cell(1, N_i);
cell_ts_1 = cell(1, N_i);
cell_Fc_1 = cell(1, N_i);
cell_Tc_1 = cell(1, N_i);

Temp_1_saved = []; %convert cell array to single vector
Conc_1_saved = [];
noise_1_saved = [];
Flow_1_saved = [];
Temp_CW_1_saved = [];
sum_N1 = 0;
time_1_saved = [];

parpool('local',12);

parfor k = 1:N_i %for loop for each crossing point at lambda_1
    
    Ts_init = Ts_saved(k); %initial value
    Cs_init = Cs_saved(k); %initial value
    Fc_init = Fc_saved(k);
    Tc_init = Tc_saved(k);
    ts_init = ts_saved(k);
    N = 0;
    Ts_1_saved = [];
    Cs_1_saved = [];
    Fc_1_saved = [];
    Tc_1_saved = [];
    dW_1_saved = [];
    ts_1_saved = [];

for j = 1:k_i
        
        Cs  = [];
        Ts = [];
        ts = [];
        Fc = [];
        Tc = [];

        Cs(1) = Cs_init;
        Ts(1) = Ts_init;
        Fc(1) = Fc_init;
        Tc(1) = Tc_init;
        ts(1) = ts_init;
        
    for i = 2:traj_length
            
    C = Cs(i-1) ;
    T = Ts(i-1) ;   
    Tc_s = Tc(i-1);
    F = (Fc_0 + Kc*(T - Tsp));
    dW = random(pd);
    noise = dW ;
  
    %% CW constraints should come here!
    if F < Fc_lb
        F = Fc_lb;
    end

    %% Create functions at current time-step
    a = A*(2-C) - ko*C*exp(-E/R/T);
    a_T = A*(300-T) + alpha*(Tc_s-T) + beta*C*exp(-E/R/T); % RHS of energy balance equation to iterate
    a_Tc = F*(300-Tc_s)/10 - alpha*(Tc_s - T);

    %% Update for next time-step
    Cs(i) = C + a*h + b*dW*sqrt(h); % Updating concentration for next time step
    
    Ts(i) = T + a_T*h  ; % Updating temperature for next time step

    Tc(i) = Tc_s + a_Tc*h;
    
    ts(i) = ts(i-1) + h;

    Fc(i) = F;

    %% Condition for crossing
     
        if  Ts(i)< Ts(i-1) && Ts(i)<lambda_i && Ts(i-1)>lambda_i %&& Cs_1(i) > solving_for_conc_as_func_of_T(tau, 510)- 1.2*sqrt(var) && Cs_1(i) < solving_for_conc_as_func_of_T(tau, 510)+ 1.2*sqrt(var)% If The trajectory crosses over the first order parameter of 580 K
         
         N = N + 1;
         Cs_1_saved = [Cs_1_saved, Cs(i)];
         Ts_1_saved = [Ts_1_saved, Ts(i)];
         dW_1_saved = [dW_1_saved, noise];
         ts_1_saved = [ts_1_saved, ts(i)];
         Fc_1_saved = [Fc_1_saved, Fc(i)];
         Tc_1_saved = [Tc_1_saved, Tc(i)];
         
         break

  
        end
       
    end
             
end

    cell_Ts_1{k} = Ts_1_saved;
   cell_Cs_1{k} = Cs_1_saved;
   cell_noise_1{k} = dW_1_saved;
    cell_N1{k} = N;
    cell_ts_1{k} = ts_1_saved;
    cell_Fc_1{k} = Fc_1_saved;
    cell_Tc_1{k} = Tc_1_saved;

    Temp_1_saved = [Temp_1_saved, Ts_1_saved];
    Conc_1_saved = [Conc_1_saved, Cs_1_saved];
    sum_N1 = sum_N1 + N;
    noise_1_saved = [noise_1_saved, dW_1_saved];  
    time_1_saved = [time_1_saved, ts_1_saved];
    Flow_1_saved = [Flow_1_saved, Fc_1_saved];
    Temp_CW_1_saved = [Temp_CW_1_saved, Tc_1_saved];

end

delete(gcp('nocreate'));

end
%%