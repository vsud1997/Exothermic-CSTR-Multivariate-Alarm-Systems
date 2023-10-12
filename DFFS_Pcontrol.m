function [overall_probability_DFFS, overall_rate_DFFS] = DFFS_Pcontrol(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Fc_ub,N_0,Rate_of_transition_1, Ts_saved,Cs_saved,Fc_saved,Tc_saved,ts_saved)

A = 1/tau; % 1/ResidenceTime =  1/0.5 = 2
b = 1/tau;
ko = 17.038; % Rate constant
E = 1.50e+04; % Actication energy
R = 8.314; % Gas Constant
alpha = 0.075; 
beta = 9370.9;
pd = makedist('Normal',0,sqrt(var)); 

h = 0.01;
N_interface = 8;
transition_time = 16; %minutes
time_traj = 3*transition_time/N_interface;
N_traj = N_0+50;
traj_length = time_traj/h + 1;
b = 1/tau;

lambda_1 = 620;
lambda_2 = 580;
lambda_3 = 548;
lambda_4 = 524;
lambda_5 = 492;
lambda_6 = 455;
lambda_7 = 425;
lambda_8 = 400;

tBounds2 = [0  time_traj];
traj_length = (tBounds2(2) - tBounds2(1))/h;

%% Compute statistics for transition from lambda_0/lambda_A to lambda_1

[Cs_1_saved, Ts_1_saved, Tc_1_saved, Fc_1_saved, N1, probability_N1] = DFFS_Pcontrol_internal_loop(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Fc_ub, lambda_1, N_0, N_traj, traj_length, Ts_saved,Cs_saved,Fc_saved,Tc_saved);

%% Compute statistics for transition from lambda_1 to lambda_2

[Cs_2_saved, Ts_2_saved, Tc_2_saved, Fc_2_saved, N2, probability_N2] = DFFS_Pcontrol_internal_loop(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Fc_ub, lambda_2, N1, N_traj, traj_length, Ts_1_saved,Cs_1_saved,Fc_1_saved,Tc_1_saved);

%% Compute statistics for transition from lambda_2 to lambda_3

[Cs_3_saved, Ts_3_saved, Tc_3_saved, Fc_3_saved, N3, probability_N3] = DFFS_Pcontrol_internal_loop(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Fc_ub, lambda_3, N2, N_traj, traj_length, Ts_2_saved,Cs_2_saved,Fc_2_saved,Tc_2_saved);

%% Compute statistics for transition from lambda_3 to lambda_4

[Cs_4_saved, Ts_4_saved, Tc_4_saved, Fc_4_saved, N4, probability_N4] = DFFS_Pcontrol_internal_loop(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Fc_ub, lambda_4, N3, N_traj, traj_length, Ts_3_saved,Cs_3_saved,Fc_3_saved,Tc_3_saved);

%% Compute statistics for transition from lambda_4 to lambda_5

[Cs_5_saved, Ts_5_saved, Tc_5_saved, Fc_5_saved, N5, probability_N5] = DFFS_Pcontrol_internal_loop(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Fc_ub, lambda_5, N4, N_traj, traj_length, Ts_4_saved,Cs_4_saved,Fc_4_saved,Tc_4_saved);
%% Compute statistics for transition from lambda_5 to lambda_6

[Cs_6_saved, Ts_6_saved, Tc_6_saved, Fc_6_saved, N6, probability_N6] = DFFS_Pcontrol_internal_loop(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Fc_ub, lambda_6, N5, N_traj, traj_length, Ts_5_saved,Cs_5_saved,Fc_5_saved,Tc_5_saved);

%% Compute statistics for transition from lambda_6 to lambda_7

[Cs_7_saved, Ts_7_saved, Tc_7_saved, Fc_7_saved, N7, probability_N7] = DFFS_Pcontrol_internal_loop(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Fc_ub, lambda_7, N6, N_traj, traj_length, Ts_6_saved,Cs_6_saved,Fc_6_saved,Tc_6_saved);

%% Compute statistics for transition from lambda_7 to lambda_7/lambda_B

[Cs_8_saved, Ts_8_saved, Tc_8_saved, Fc_8_saved, N8, probability_N8] = DFFS_Pcontrol_internal_loop(tau, var, Kc, Tsp, Tc_0, Fc_0, Fc_lb, Fc_ub, lambda_8, N7, N_traj, traj_length, Ts_7_saved,Cs_7_saved,Fc_7_saved,Tc_7_saved);

overall_probability_DFFS = probability_N1*probability_N2*probability_N3*probability_N4*probability_N5*probability_N6*probability_N7*probability_N8;

overall_rate_DFFS = Rate_of_transition_1*overall_probability_DFFS; % min-1



end