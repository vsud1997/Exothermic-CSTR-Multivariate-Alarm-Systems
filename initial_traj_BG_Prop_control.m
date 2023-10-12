function [Init_traj_output] = initial_traj_BG_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0, tBounds, Ntraj, lambda_A)

A = 1/tau; % 1/ResidenceTime =  1/0.5 = 2
ko = 17.038; % Rate constant
E = 1.50e+04; % Actication energy
R = 8.314; % Gas Constant
alpha = 0.075; 
beta = 9370.9;
% Step size and other details
h = 0.01;  %step size
N = (tBounds(2)-tBounds(1))/h; % Number of steps
%b = 2;
b = 1/tau;
%ts = linspace(tBounds(1), tBounds(2), N); % From t0-->t1 with N points
pd = makedist('Normal',0,sqrt(var)); % Probability distribution of mean 0 and variance dt
Selection_Array_1 = []; %Numbering from 1 till length of total crossings

time_instants_basin_A = [];

C_init = 1.25; %kmol/m3
T_init = 660; %K
Tc_init = Tc_0;

N_0 = 0; %Number of crossing points for initial trajectory

Ts_saved = [];
noise_saved = [];
ts_saved = [];
Cs_saved = [];
Tc_saved = [];
Fc_saved = [];

ts_prior = [];
Cs_prior = [];
Ts_prior = [];
Tc_prior = [];
Fc_prior = [];


flag = 0;


for j = 1:Ntraj

    Cs  = zeros(1,N); %Defining array to store computed Concentrations
    Ts  = zeros(1,N); %Defining array to store computed temperatures
    Tc_s = zeros(1,N);
    Fc_s = zeros(1,N);
    ts = zeros(1,N);

    Cs(1) = C_init; % Initial value of concentration
    Ts(1) = T_init; % Initial value of temperature
    Tc_s(1) = Tc_init;
    Fc_s(1) = Fc_0;
    
    
    for i = 2:N

    C = Cs(i-1) ;
    T = Ts(i-1) ;
    Tc = Tc_s(i-1) ;
    F = (Fc_0 + Kc*(T - Tsp));
    t = ts(i-1) ;
    dW = random(pd);
    noise = dW ;

     %% Create functions at current time-step
    a = A*(2-C) - ko*C*exp(-E/R/T);
    a_T = A*(300-T) + alpha*(Tc-T) + beta*C*exp(-E/R/T); % RHS of energy balance equation to iterate
    a_Tc = F*(Tc_0-Tc)/10 - alpha*(Tc - T);

    %% Update for next time-step
    
    Cs(i) = C + a*h + b*dW*sqrt(h); % Updating concentration for next time step
    
    Ts(i) = T + a_T*h  ; % Updating temperature for next time step
    
    Tc_s(i) = Tc + a_Tc*h ;
    
    Fc_s(i) = F;
        
    ts(i) = t + h ;
    
    
     if Ts(i)>500
         
         time_instants_basin_A = [time_instants_basin_A, ts(i)];
         
     end
    
     if  Ts(i)<Ts(i-1) && Ts(i)<=lambda_A && Ts(i-1)>lambda_A %&& Cs(i) > solving_for_conc_as_func_of_T(tau, 650)- 2.5*sqrt(var) && Cs(i) < solving_for_conc_as_func_of_T(tau, 650)+ 2.5*sqrt(var)  % If The trajectory crosses over the first order parameter
         
         N_0 = N_0 + 1;
         Ts_saved = [Ts_saved, Ts(i)] ;% Storing the particular temperature values;
         Cs_saved = [Cs_saved, Cs(i)] ; %Storing the particular concentration values;
         noise_saved = [noise_saved, noise];
         ts_saved = [ts_saved, ts(i)];
         Tc_saved = [Tc_saved, Tc_s(i)];
         Fc_saved = [Fc_saved, Fc_s(i)];


     end
         
     if N_0 >= 100
         
         flag = 1;
         break
         
     end
    end
    
    %m = find(Ts>500);
    %time_instants_basin_A = [time_instants_basin_A, ts(m)];
    
   ts_prior_cell{j} = ts_prior;
   Ts_prior_cell{j} = Ts_prior;
   Cs_prior_cell{j} = Cs_prior;
   Fc_prior_cell{j} = Fc_prior;
   Tc_prior_cell{j} = Tc_prior;

    
    if flag==1
        break
    end


end

sum = 0;

time_instants_basin_A = time_instants_basin_A/h;

for j = 2: length(time_instants_basin_A)
    
    if  (time_instants_basin_A(j) - time_instants_basin_A(j-1) ) == 1
        
        sum = sum + 1;
    
    end
        
end

time_basin_A = sum*0.01;

Initial_Rate_of_transition_1 = N_0/time_basin_A;

time_instants_basin_A = 0.01*time_instants_basin_A;

n_output = 8; %number of output variables 
Init_traj_output = cell(n_output, 1);

Init_traj_output{1} = N_0;
Init_traj_output{2} = Ts_saved;
Init_traj_output{3} = Cs_saved;
Init_traj_output{4} = noise_saved;
Init_traj_output{5} = Tc_saved;
Init_traj_output{6} = Fc_saved;
Init_traj_output{7} = ts_saved;
Init_traj_output{8} = Initial_Rate_of_transition_1;


end