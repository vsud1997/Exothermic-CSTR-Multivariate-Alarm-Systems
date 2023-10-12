function [N_L, N_ESD, N_L_to_ESD, N_basin_B, time_L_ESD, time_ESD_basin_B, Avg_N_B ] = dynamic_risk_analysis_unpostulated(N_sim, tau_init)

V = 10;
ko = 17.038; % Rate constant
E = 1.50e+04; % Activation energy
R = 8.314; % Gas Constant
alpha = 0.075; 
beta = 9370.9;
Fc0 = 50; %m3/min % Flow-rate of CW at steady state
Tsp = 800; %k
Tc0 = 300; %K Initial temperature of CW
Kc = 0.02;% Controller gain
% Step size and other details
tBounds = [0 70000];        % The bounds of t
var = 0.02; % Variance
h = 0.01;  %step size

N = (tBounds(2)-tBounds(1))/h + 1; % Number of steps

C_init = 1.2; %kmol/m3
T_init = 700; %K
Tc_init = Tc0; %K
Fc_init = Fc0;
pd = makedist('Normal',0,sqrt(var)); % Probability distribution of mean 0 and variance dt

ts    = linspace(tBounds(1), tBounds(2), N); % From t0-->t1 with N points


cell_Ts = cell(N_sim,1);
cell_Fc = cell(N_sim,1);
cell_Tc = cell(N_sim,1);
cell_Cs = cell(N_sim,1);


Ts_L = 501;
Fc_L = 44.05;
Tc_L = 303.4;

%Ts_LL = 509.5;
%Fc_LL = 44.21;
%Tc_LL = 303.52;

%Ts_LLL = 487.5;
%Fc_LLL = 43.78;
%Tc_LLL= 303.19;

Ts_ESD = 453.1;
Fc_ESD = 43.1;
Tc_ESD = 302.6572;


false_positive_ESD = 0;
true_positive_ESD = 0;


cell_time_from_ESD_to_basin_B = cell(1,N_sim);

cell_Simulation_number_false_positive = cell(N_sim,1);
cell_Simulation_number_true_positive = cell(N_sim,1);
cell_time_ESD_per_sim = cell(N_sim,1);
cell_timestamp_basin_B = cell(N_sim,1);

cell_tau = cell(N_sim,1);
cell_N_A_0 = cell(N_sim,1);
cell_N_A = cell(N_sim,1);
cell_N_B = cell(N_sim,1);

ESD_first_time = 0;
basin_B_first_time = 0;

N_ESD = 0;
N_L = 0;

N_basin_B = 0;

cell_fall_time_instants = cell(1,N_sim);
cell_rise_time_instants = cell(1, N_sim);

 cell_time_spent_in_B = cell(1,N_sim);

cell_time_basin_B = cell(1,N_sim);

L_tag = [];
ESD_tag = [];
ESD_source_tag = [];
L_tag_counter = 0;
ESD_tag_counter = 0;


% define limits for uniformly distributed response times and response
% actions
op_timestep_lower = 23; %response time of 15 seconds
op_timestep_upper = 155; %response time of 60 seconds
op_action_lower = 0.5346; %least aggressive response,i.e., no action taken
op_action_upper = 0.5690; %most aggressive response

cell_timestamps_op_action = cell(1, N_sim);

% define arrays for operator response times and actions

op_response_time = [];
op_action = [];


for j = 1:N_sim

tau = tau_init;
time_spent_in_basin_B = [];
Simulation_number_false_positive = [];
Simulation_number_true_positive = [];

time_ESD_per_sim = [];    



time_L_per_sim = [];

counter_basin_B = 0;




Cs  = []; %Defining array to store computed Concentrations
Ts  = []; %Defining array to store computed temperatures
Tc_s = [];%Defining array to store computed  CW temperatures
Fc_s = []; %Defining array to store computed  CW flow-rate
ts = [];


Cs(1) = C_init; % Initial value of concentration
Ts(1) = T_init; % Initial value of temperature
Tc_s(1) = Tc_init;% Initial value of CW temperature
Fc_s(1) = Fc_init ;
ts(1) = 0;

alarm_k_L = 0; %counter for low alarm 1 
alarm_k_L = 0;
alarm_k_LLL = 0;
alarm_k_ESD = 0;

timestamps_L_on = [];
timestamps_LL_on = [];
timestamps_L_off = [];
timestamps_LL_off = [];

time_from_L_to_LL = [];

L_alarm_tag = [];
LL_alarm_tag = [];
LL_alarm_source_tag = [];


time_from_LL_to_LLL  = [];
time_from_L_to_ESD = [];
time_from_ESD_to_basin_B = [];
time_ESD_basin_B = 0;
time_spent_in_basin_B = [];

N_A = [];
N_A_0 = [];
N_B = [];
tau_per_sim = [];

N_A_0(1) = V*2/tau_init;
N_A(1) = V*C_init/tau_init;
N_B(1) =  V*(2 - C_init)/tau_init;
tau_per_sim(1) = tau_init;
% define array for operator response timestamps
op_response_timestamp = [];
t_alarm = 0;

for i = 2:N
    
    % Stochastic DE given as : dX = a(X)dt + b(X)dW where X: Stochastic
    % Process and dW: Noise
    A = 1/tau;
    b = 1/tau;
    C = Cs(i-1) ;
    T = Ts(i-1) ;
    Tc = Tc_s(i-1);
    a = A*(2-C) - ko*C*exp(-E/R/T);
    
    
    
    dW = random(pd);
        
  
   
    a_T = A*(300-T) + alpha*(Tc-T) + beta*C*exp(-E/R/T); % RHS of energy balance equation to iterate
    
    a_Tc = (Fc0 + Kc*(T - Tsp))*(300-Tc)/10 - alpha*(Tc - T);
    
    Fc_s(i) = Fc0 + Kc*(T - Tsp);
    
    a_error = T - Tsp ;
    
    Cs(i) = C + a*h + b*dW*sqrt(h)  ; % Updating concentration for next time step
    
    Ts(i) = T + a_T*h  ; % Updating temperature for next time step
    
    Tc_s(i) = Tc + a_Tc*h ;

    ts(i) = ts(i-1) + h;

    %%% If-else test for start time and end time of alarms###
    %%% Part 1: LOW level alarms ####

    if Ts(i)<=Ts_L && Fc_s(i)<=Fc_L && Tc_s(i)<=Tc_L && alarm_k_L == 0

        alarm_k_L = 1;
        time_L = ts(i);
        time_L_per_sim = [time_L_per_sim, time_L];
        N_L = N_L + 1;
        L_tag_counter = L_tag_counter + 1;
        L_tag = [L_tag, L_tag_counter];

        op_response_time_var = randi([op_timestep_lower, op_timestep_upper] ,1,1); %has to be an integer
        op_action_var = op_action_lower + (op_action_upper - op_action_lower)*rand(1,1);
        t_alarm = i;

        % save the operator response time and action
        op_response_time = [op_response_time, op_response_time_var];
        op_action = [op_action, op_action_var];

        %tau = op_action_var;
        
    elseif Ts(i)>Ts_L && Fc_s(i)>Fc_L && Tc_s(i)>Tc_L && alarm_k_L == 1

        alarm_k_L = 0;
        tau = tau_init;

    end

    if t_alarm ~=0 && i == t_alarm + op_response_time_var && alarm_k_ESD == 0 % time time when the action will be taken

        tau = op_action_var;
        op_response_timestamp = [op_response_timestamp, ts(i)];
        

    end


    if Ts(i)<=Ts_ESD && Fc_s(i)<=Fc_ESD && Tc_s(i)<=Tc_ESD && alarm_k_ESD == 0

        alarm_k_ESD = 1;
        time_ESD_basin_B = ts(i);
        time_ESD_per_sim = [time_ESD_per_sim, time_ESD_basin_B];
        N_ESD = N_ESD + 1;
        time_from_L_to_ESD = [time_from_L_to_ESD, time_ESD_basin_B - time_L];
        ESD_tag_counter = ESD_tag_counter + 1;
        ESD_tag = [ESD_tag, ESD_tag_counter];
        ESD_source_tag = [ESD_source_tag, L_tag_counter];

        tau = 0.62;

        elseif Ts(i)>Ts_ESD && Fc_s(i)>Fc_ESD && Tc_s(i)>Tc_ESD && alarm_k_ESD == 1

        alarm_k_ESD = 0;
        tau = tau_init;
        
    end

    if Ts(i)<=400 && counter_basin_B == 0 
            
            timestamp_basin_B = ts(i);
            time_from_ESD_to_basin_B = [time_from_ESD_to_basin_B, timestamp_basin_B - time_ESD_basin_B];
            counter_basin_B = 1;
            basin_B_first_time = 1;
            N_basin_B = N_basin_B + 1;

    elseif Ts(i)>440 && counter_basin_B == 1

        counter_basin_B = 0;
        time_spent_in_basin_B = [time_spent_in_basin_B, ts(i) - timestamp_basin_B];

    end

        N_A_0(i) = 2*V/tau;
        N_A(i) = V*Cs(i)/tau;
        N_B(i) = V*(2 - Cs(i))/tau;
        tau_per_sim(i) = tau;
end

if isempty(time_from_ESD_to_basin_B) == 0

    true_positive_ESD = true_positive_ESD + 1;
    Simulation_number_true_positive = [Simulation_number_true_positive, j];

end

[f_fall,lt_fall,ut_fall,ll_fall,ul_fall] = falltime(Ts,ts);
[r_rise,lt_rise,ut_rise,ll_rise,ul_rise] = risetime(Ts,ts);

if length(lt_fall)<=5

    cell_fall_time_instants{j} = lt_fall;
    cell_rise_time_instants{j} = lt_rise;

end

cell_time_ESD_per_sim{j} = time_ESD_per_sim;
%cell_time_LLL_per_sim{j} = time_LLL_per_sim;
cell_time_L_per_sim{j} = time_L_per_sim;


cell_time_from_ESD_to_basin_B{j} = time_from_ESD_to_basin_B(find(time_from_ESD_to_basin_B<=10));
cell_time_from_L_to_ESD{j} = time_from_L_to_ESD(find(time_from_L_to_ESD<15));

cell_Ts{j} = Ts;
cell_Fc{j} = Fc_s;
cell_Tc{j} = Tc_s;
cell_Cs{j} = Cs;

cell_Simulation_number_false_positive{j} = Simulation_number_false_positive;
cell_Simulation_number_true_positive{j} = Simulation_number_true_positive;

cell_time_basin_B{j} = time_spent_in_basin_B;

k2 = find(Cs>=0 & Cs<=2); %consider only the realistic concentrations

cell_N_A{j} = mean(N_A(k2));
cell_N_A_0{j} = mean(N_A_0(k2));
cell_N_B{j} = mean(N_B(k2));
cell_tau{j} = tau_per_sim(k2);

cell_timestamps_op_action{j} = op_response_timestamp;

end

for k = 1:N_sim
    if isempty(cell_fall_time_instants{k})== 0
        rise_time = cell_rise_time_instants{k};
        fall_time = cell_fall_time_instants{k};
      if length(fall_time) == length(rise_time)
        cell_time_spent_in_B{k} = rise_time - fall_time;
      else
        cell_time_spent_in_B{k} = rise_time - fall_time(1:end-1);
      end
    end
 end

unique_ESD_source_tag = unique(ESD_source_tag);
N_L_to_ESD = length(unique_ESD_source_tag);

time_L_ESD = cell2mat(cell_time_from_L_to_ESD);
time_ESD_basin_B = cell2mat(cell_time_from_ESD_to_basin_B);

Avg_N_B = mean(cell2mat(cell_N_B));

end
