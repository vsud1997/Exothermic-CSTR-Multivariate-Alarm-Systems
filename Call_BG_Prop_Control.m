clc
clearvars
tic
var = 0.02;
tBounds = [0 30000];
Ntraj = 10;

tau = 0.53;
Kc = 0.06;

Fc_0 = 50;
Fc_ub = 70;
Fc_lb = 30; 
Tsp = 800;
Tc_0 = 300;

k0 = 70;
k1 = 60;
k2 = 50;
k3 = 40;
k4 = 40;
k5 = 40;
k6 = 40;
k7 = 40;



total = k0*k1*k2*k3*k4*k5*k6*k7;
lambda_A = 650;

%% Initial Rate analyses as function of Kc and tau
%Initial_Rate_array = [];
%parpool('local',10);
%parfor j = 1:25
[Init_traj_output] = initial_traj_BG_Prop_control(tau, var, Kc, Tsp, Tc_0, Fc_0, tBounds, Ntraj,lambda_A);
N_0 = Init_traj_output{1};
Ts_saved = Init_traj_output{2};
Cs_saved = Init_traj_output{3};
noise_saved = Init_traj_output{4};
Tc_saved = Init_traj_output{5};
Fc_saved = Init_traj_output{6};
ts_saved = Init_traj_output{7};
Initial_Rate_of_transition_1 = Init_traj_output{8};
%Initial_Rate_array = [Initial_Rate_array, Init_traj_output{8}];
%end

%mean(Initial_Rate_array)
toc
%delete(gcp('nocreate'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%Store initial configurations corresponding to respective committor %probabilities

%%%Variables at lambda_0
Cs_0_init = Cs_new_guess;
Ts_0_init = Ts_new_guess;
noise_0_init = noise_new_guess;
ts_0_init = ts_new_guess;
Fc_0_init = Fc_new_guess;
Tc_0_init = Tc_new_guess;

%%%% Variables at lambda_1
Cs_1_init = Cs_1_saved;
Ts_1_init = Ts_1_saved;
noise_1_init = dW_1_saved;
ts_1_init = ts_1_saved;
Fc_1_init = Fc_1_saved;
Tc_1_init = Tc_1_saved;

%%%% Variables at lambda_2
Cs_2_init = cell2mat(cell_Cs_1);
Ts_2_init = cell2mat(cell_Ts_1);
noise_2_init = cell2mat(cell_noise_1);
%ts_2_init = cell2mat(cell_ts_1);
ts_2_init = cell2mat(cell_ts_1);     
Fc_2_init = cell2mat(cell_Fc_1);
Tc_2_init = cell2mat(cell_Tc_1);

%%%% Variables at lambda_3
Cs_3_init = cell2mat(cell_Cs_2);
Ts_3_init = cell2mat(cell_Ts_2);
noise_3_init = cell2mat(cell_noise_2);
%ts_3_init = cell2mat(cell_ts_2);
ts_3_init = cell2mat(cell_ts_2);
Fc_3_init = cell2mat(cell_Fc_2);
Tc_3_init = cell2mat(cell_Tc_2);

%%%% Variables at lambda_4
Cs_4_init = cell2mat(cell_Cs_3);
Ts_4_init = cell2mat(cell_Ts_3);
noise_4_init = cell2mat(cell_noise_3);
%ts_4_init = cell2mat(cell_ts_3);
ts_4_init = cell2mat(cell_ts_3);
Fc_4_init = cell2mat(cell_Fc_3);
Tc_4_init = cell2mat(cell_Tc_3);


%%%% Variables at lambda_5
Cs_5_init = cell2mat(cell_Cs_4);
Ts_5_init = cell2mat(cell_Ts_4);
noise_5_init = cell2mat(cell_noise_4);
%ts_5_init = cell2mat(cell_ts_4);
ts_5_init = cell2mat(cell_ts_4);
Fc_5_init = cell2mat(cell_Fc_4);
Tc_5_init = cell2mat(cell_Tc_4);


%%%% Variables at lambda_6
Cs_6_init = cell2mat(cell_Cs_5);
Ts_6_init = cell2mat(cell_Ts_5);
noise_6_init = cell2mat(cell_noise_5);
%ts_6_init = cell2mat(cell_ts_5);
ts_6_init = cell2mat(cell_ts_5);
Fc_6_init = cell2mat(cell_Fc_5);
Tc_6_init = cell2mat(cell_Tc_5);

%%%% Variables at lambda_7
Cs_7_init = cell2mat(cell_Cs_6);
Ts_7_init= cell2mat(cell_Ts_6);
noise_7_init = cell2mat(cell_noise_6);
%ts_6_init = cell2mat(cell_ts_5);
ts_7_init = cell2mat(cell_ts_6);
Fc_7_init = cell2mat(cell_Fc_6);
Tc_7_init = cell2mat(cell_Tc_6);

toc


%%%% Compute Max and Min CW flow-rates

min_Fc_array = []; 
min_Fc_array = [min_Fc_array, [Fc_0_init,min(Fc_1_init),min(Fc_2_init),min(Fc_3_init),min(Fc_4_init),min(Fc_5_init),min(Fc_6_init),min(Fc_7_init)]];

max_fc_array = []; 
max_fc_array = [max_fc_array, [Fc_0_init,max(Fc_1_init),max(Fc_2_init),max(Fc_3_init),max(Fc_4_init),max(Fc_5_init),max(Fc_6_init),max(Fc_7_init)]];
    
%%% Call committor prob function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%[p_B_6, p_B_5, p_B_4, p_B_3, p_B_2, p_B_1, p_B_0] = p_B_committor(cell_N1, cell_N2, cell_N3,cell_N4,cell_N5,cell_N6, k0, k1,k2,k3,k4,k5,k6);

[p_B_7, p_B_6, p_B_5, p_B_4, p_B_3, p_B_2, p_B_1, p_B_0] = p_B_committor(cell_N1, cell_N2, cell_N3,cell_N4,cell_N5,cell_N6, cell_N7, k0, k1,k2,k3,k4,k5,k6, k7);

toc

%%% Committor Probabilities at lambda_7 = 1

%p_B_7 = ones(1, length(Ts_7_init));

committor_prob = [];
committor_prob = [committor_prob, p_B_0];
committor_prob = [committor_prob, transpose(p_B_1)];
committor_prob = [committor_prob, transpose(p_B_2)];
committor_prob = [committor_prob, transpose(p_B_3)];
committor_prob = [committor_prob, transpose(p_B_4)];
committor_prob = [committor_prob, transpose(p_B_5)];
committor_prob = [committor_prob, transpose(p_B_6)];
committor_prob = [committor_prob, p_B_7];
toc

tic
figure(1)
%h = histogram(committor_prob_2_sims)
%hist_current_sim = histogram(Tall_p_B)
histogram(p_B_7, 10)
title(" Histogram of committor probabilities")
xlabel('p_B')
ylabel('Occurrences')
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compress committer data to more
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% meaningful points %%%%

k_compress_pB_3 = find(p_B_3 > mean(p_B_3)- 0.5*std(p_B_3) & p_B_3 < mean(p_B_3)+ 0.5*std(p_B_3));
k_compress_pB_4 = find(p_B_4 > mean(p_B_4)- 0.5*std(p_B_4) & p_B_4 < mean(p_B_4)+ 0.5*std(p_B_4));
k_compress_pB_5 = find(p_B_5 > mean(p_B_5)- 0.5*std(p_B_5) & p_B_5 < mean(p_B_5)+ 0.5*std(p_B_5));
k_compress_pB_6 = find(p_B_6 > mean(p_B_6)- 0.5*std(p_B_6) & p_B_6 < mean(p_B_6)+ 0.5*std(p_B_6));
k_compress_pB_7 = find(p_B_7 > mean(p_B_7) - 0.1*std(p_B_7) &  p_B_7 < mean(p_B_7) + 0.1*std(p_B_7));
%k_compress_pB_7 = randperm(length(Ts_7_init), floor(0.1*length(p_B_6)));

p_B_3_compressed = p_B_3(k_compress_pB_3);
p_B_4_compressed = p_B_4(k_compress_pB_4);
p_B_5_compressed = p_B_5(k_compress_pB_5);
p_B_6_compressed = p_B_6(k_compress_pB_6);
p_B_7_compressed = p_B_7(k_compress_pB_7);

%%%Compress all variables at lambda_7

committer_prob_final = [];
committer_prob_final = [committer_prob_final, p_B_0];
committer_prob_final = [committer_prob_final, transpose(p_B_1)];
committer_prob_final = [committer_prob_final, transpose(p_B_2)];
committer_prob_final = [committer_prob_final, transpose(p_B_3_compressed)];
committer_prob_final = [committer_prob_final, transpose(p_B_4_compressed)];
committer_prob_final = [committer_prob_final, transpose(p_B_5_compressed)];
committer_prob_final = [committer_prob_final, transpose(p_B_6_compressed)];
committer_prob_final = [committer_prob_final, p_B_7_compressed];

%%%% Do the same for the process variables %%%%%%%%

Cs_3_compressed = Cs_3_init(k_compress_pB_3);
Ts_3_compressed = Ts_3_init(k_compress_pB_3);
noise_3_compressed = noise_3_init(k_compress_pB_3);
ts_3_compressed = ts_3_init(k_compress_pB_3);
Fc_3_compressed = Fc_3_init(k_compress_pB_3);
Tc_3_compressed = Tc_3_init(k_compress_pB_3);

Cs_4_compressed = Cs_4_init(k_compress_pB_4);
Ts_4_compressed = Ts_4_init(k_compress_pB_4);
noise_4_compressed = noise_4_init(k_compress_pB_4);
ts_4_compressed = ts_4_init(k_compress_pB_4);
Fc_4_compressed = Fc_4_init(k_compress_pB_4);
Tc_4_compressed = Tc_4_init(k_compress_pB_4);

Cs_5_compressed = Cs_5_init(k_compress_pB_5);
Ts_5_compressed = Ts_5_init(k_compress_pB_5);
noise_5_compressed = noise_5_init(k_compress_pB_5);
ts_5_compressed = ts_5_init(k_compress_pB_5);
Fc_5_compressed = Fc_5_init(k_compress_pB_5);
Tc_5_compressed = Tc_5_init(k_compress_pB_5);

Cs_6_compressed = Cs_6_init(k_compress_pB_6);
Ts_6_compressed = Ts_6_init(k_compress_pB_6);
noise_6_compressed = noise_6_init(k_compress_pB_6);
ts_6_compressed = ts_6_init(k_compress_pB_6);
Fc_6_compressed = Fc_6_init(k_compress_pB_6);
Tc_6_compressed = Tc_6_init(k_compress_pB_6);

Cs_7_compressed = Cs_7_init(k_compress_pB_7);
Ts_7_compressed = Ts_7_init(k_compress_pB_7);
noise_7_compressed = noise_7_init(k_compress_pB_7);
ts_7_compressed = ts_7_init(k_compress_pB_7);
Fc_7_compressed = Fc_7_init(k_compress_pB_7);
Tc_7_compressed = Tc_7_init(k_compress_pB_7);

Cs_overall = [];
Ts_overall = [];
noise_overall = [];
time_overall = [];
Fc_overall = [];
Tc_overall = [];

Cs_overall = [Cs_overall, Cs_0_init, Cs_1_init, Cs_2_init, Cs_3_compressed, Cs_4_compressed, Cs_5_compressed, Cs_6_compressed, Cs_7_compressed]; 
Ts_overall = [Ts_overall, Ts_0_init, Ts_1_init, Ts_2_init, Ts_3_compressed, Ts_4_compressed, Ts_5_compressed, Ts_6_compressed, Ts_7_compressed]; 
noise_overall = [noise_overall,  noise_0_init, noise_1_init, noise_2_init, noise_3_compressed,noise_4_compressed, noise_5_compressed, noise_6_compressed, noise_7_compressed];
time_overall = [time_overall, ts_0_init, ts_1_init, ts_2_init, ts_3_init, ts_4_compressed, ts_5_compressed, ts_6_compressed, ts_7_compressed];
Fc_overall = [Fc_overall, Fc_0_init, Fc_1_init, Fc_2_init, Fc_3_compressed, Fc_4_compressed, Fc_5_compressed, Fc_6_compressed, Fc_7_compressed];
Tc_overall = [Tc_overall, Tc_0_init, Tc_1_init, Tc_2_init, Tc_3_compressed, Tc_4_compressed, Tc_5_compressed, Tc_6_compressed, Tc_7_compressed];



%%%% Probability distribution histogram for T = 620 K %%%%%% 
figure(2)
histogram(p_B_7_compressed,35)
title(" Histogram of committor probabilities for T = 620 K")
xlabel('p_B')
ylabel('Occurrences')

%%%% Probability distribution histogram for T = 580 K %%%%%% 
figure(3)
histogram(p_B_2)
title(" Histogram of committor probabilities for T = 580 K")
xlabel('p_B')
ylabel('Occurrences')


%%%% Probability distribution histogram for T = 548 K %%%%%% 
figure(4)
histogram(p_B_3_compressed)
title(" Histogram of committor probabilities for T = 548 K")
xlabel('p_B')
ylabel('Occurrences')

%%%% Probability distribution histogram for T = 524 K %%%%%% 
figure(5)
histogram(p_B_4_compressed,10)
title(" Histogram of committor probabilities for T = 524 K")
xlabel('p_B')
ylabel('Occurrences')

%%%% Probability distribution histogram for T = 492 K %%%%%% 
figure(6)
histogram(p_B_5_compressed,10)
title(" Histogram of committor probabilities for T = 492 K")
xlabel('p_B')
ylabel('Occurrences')

%%%% Probability distribution histogram for T = 455 K %%%%%% 
figure(7)
histogram(p_B_7_compressed)
title(" Histogram of committor probabilities for T = 455 K")
xlabel('p_B')
ylabel('Occurrences')


%%%%%%%%%%%%%%% Plot comiittor as a function of Fc%%%%%%%%%%%%%%
figure(8)
plot(Fc_overall, committer_prob_final, 'ro')
title('P_B as function of CW flow-rate, F_c')
xlabel('F_c (m^3/min)')
ylabel('Committer Probability, p_B')

%%%%%%%%%%%%%% Plot committor as a function of Tc%%%%%%%%%%
figure(9)
plot(Tc_overall, committer_prob_final, 'ko')
title('P_B as function of CW temperature, T_c')
xlabel('T_c (K)')
ylabel('Committer Probability, p_B')


%%%%%%%%%%%%%% Plot committor as a function of C%%%%%%%%%%
figure(10)
plot(Cs_overall, committer_prob_final, 'bo')
title('P_B as function of concentration')
xlabel('C_A (kmol/m^3)')
ylabel('Committer Probability, p_B')

%%%%%%%%%%%%%% Plot committor as a function of T%%%%%%%%%%
figure(11)
plot(Ts_overall, committer_prob_final, 'ko')
title('p_B as a function of temperature')
xlabel('Temp (K)')
ylabel('Committer probabilities, p_B')

%%%%%%%%%%%%%% Plot committor as a function of time %%%%%%%%%%
figure(12)
plot(time_overall,'ko')
title('p_B as a function of temperature')
xlabel('time (min)')
ylabel('Committer probabilities, p_B')



