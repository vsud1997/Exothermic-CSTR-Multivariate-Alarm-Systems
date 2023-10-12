clc
X1 = transpose(Ts_overall);
X2 = transpose(Fc_overall);
X3 = transpose(Tc_overall);
X4 = transpose(Cs_overall);
ydata = transpose(committer_prob_final);

k = find(ydata == 0) %%remove points with exact zero probabilities
ydata(k) = [];
X1(k) = []; X2(k) = []; X3(k) = []; X4(k) = [];

Xdata = [X1 X2 X3 X4];
ydata = log(ydata);

X1_rescaled = rescale(X1);
X2_rescaled = rescale(X2);
X3_rescaled = rescale(X3);
X4_rescaled = rescale(X4);

Xdata_rescaled = [X1_rescaled X2_rescaled X3_rescaled X4_rescaled];

%Split into training set and testing set
[m,n] = size(Xdata_rescaled);
P = 0.8;
idx = randperm(m);
Xdata_train = Xdata_rescaled(idx(1:round(P*m)),:);
Xdata_train_unscaled = Xdata(idx(1:round(P*m)),:);
Xdata_test = Xdata_rescaled(idx(round(P*m)+1:end),:);
Xdata_test_unscaled = Xdata(idx(round(P*m)+1:end),:);

ydata_train = ydata(idx(1:round(P*m)),:);
ydata_test = ydata(idx(round(P*m)+1:end),:);

%Define the proposed model
modelfun = @(b,x)b(1) + b(2)*x(:,1)+ b(3)*x(:,2) + b(4)*x(:,3) + b(5)*x(:,4);

b0 = [0 -4 -3 -2 -1];

tic
[b_est_lsq, resnorm, residuals] = lsqcurvefit(modelfun, b0, Xdata_train, ydata_train);
mdl = fitnlm(Xdata_train, ydata_train, modelfun, b0);
toc

%y_pred = b_est(1).*exp(b_est(2)*Xdata_rescaled(:,1)) + b_est(3).*exp(b_est(4)*Xdata_rescaled(:,2)) + b_est(5).*exp(b_est(6)*Xdata_rescaled(:,3));

b_est = mdl.Coefficients.Estimate;
y_fit = mdl.Fitted;
toc

y_mdl_train_data = modelfun(b_est_lsq, Xdata_train);
y_mdl_test_data = modelfun(b_est_lsq, Xdata_test);

%RMSPE = sqrt(mean((y_mdl_test_data - ydata_test).^2))
RMSTE_log = sqrt(mean((y_mdl_train_data - ydata_train).^2))
RMSPE_log = sqrt(mean((y_mdl_test_data - ydata_test).^2))

RMSTE_linear = sqrt(mean((exp(y_mdl_train_data) - exp(ydata_train)).^2))
RMSPE_linear = sqrt(mean((exp(y_mdl_test_data) - exp(ydata_test)).^2))

R2_log = mdl.Rsquared.Ordinary

%%%% Compute R2 Manually
SSE_array = [];
SST_array = [];
SSE_array_linear = [];
SST_array_linear = [];

tic
parpool('local',12);

mean_ydata_train = mean(ydata_train);

parfor i = 1:length(ydata_train)

    SSE_array = [SSE_array, (ydata_train(i) - y_mdl_train_data(i))^2];

    SST_array = [SST_array, (ydata_train(i) - mean_ydata_train)^2];

    SSE_array_linear = [SSE_array_linear, (exp(ydata_train(i)) - exp(y_mdl_train_data(i)))^2];

    SST_array_linear = [SST_array_linear, (exp(ydata_train(i)) - exp(mean_ydata_train))^2];



end

SSE_log = sum(SSE_array)
SST_log = sum(SST_array)

SSE_linear = sum(SSE_array_linear)
SST_linear = sum(SST_array_linear)

R2_log_manual = 1 - SSE_log/SST_log

R2_linear = 1 - SSE_linear/SST_linear

Training_error_log_model = RMSTE_linear/mean(exp(ydata_train))
Testing_error_log_model = RMSPE_linear/mean(exp(ydata_test))

toc
delete(gcp('nocreate'));

figure(1)
plot(Xdata_train_unscaled(:,1), ydata_train, 'ko')
hold on
plot(Xdata_train_unscaled(:,1), y_mdl_train_data,'ro');
%hold on
%errorbar(X1,y_fit)
hold off
legend('p_B Data', 'p_B fit', 'Location','best')
xlabel('Temperature (K)')
ylabel('ln(p_B)')
title('ln(p_B) as function of temperature, Data vs. Model Fit')

figure(3)
plot(Xdata_train(:,2), ydata_train, 'ko')
hold on
plot(Xdata_train(:,2), y_fit,'ro');
hold off
legend('p_B Data', 'p_B fit', 'Location','best')
xlabel('CW Flow-rate (m^3/min)')
ylabel('Committer Probabilty, p_B')
title('p_B as function of CW Flow-rate, Data vs. Model Fit')

figure(4)
plot(Xdata_train(:,3), ydata_train, 'ko')
hold on
plot(Xdata_train(:,3), y_fit,'ro');
hold off
legend('p_B Data', 'p_B fit', 'Location','best')
xlabel('CW temperature (K)')
ylabel('Concentration (kmol/m^3)')
title('p_B as function of CW Temperature, Data vs. Model Fit')


figure(5)
plot(Xdata_train(:,4), ydata_train, 'ko')
hold on
plot(Xdata_train(:,4), y_fit,'ro');
hold off
legend('p_B Data', 'p_B fit', 'Location','best')
xlabel('Concentration (kmol/m^3)')
ylabel('Committer Probabilty, p_B')
title('p_B as function of Concentration, Data vs. Model Fit')


corrcoef_pB_T = corrcoef(X1, ydata);
corrcoef_pB_Fc = corrcoef(X2, ydata);
corrcoef_pB_Tc = corrcoef(X3, ydata);
corrcoef_pB_C = corrcoef(X4, ydata);


%%%%%%%%%%%%%%%%
%gaussian_model = fitrgp(transpose(Cs_overall), ydata);

%%%%%%%%%%%%%%%%%%% Critical values estimation for the variables
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%

N_points = 100000;

x1_sample_unscaled = transpose(linspace(min(X1), max(X1), N_points));
x2_sample_unscaled = transpose(linspace(min(X2), max(X2), N_points));
x3_sample_unscaled = transpose(linspace(min(X3), max(X3), N_points));
x4_sample_unscaled = transpose(linspace(max(X4), min(X4), N_points));

%x1_sample_unscaled = (max(X1) - min(X1)).*rand(N_points,1) + min(X1);
%x2_sample_unscaled = (max(X2) - min(X2)).*rand(N_points,1) + min(X2);
%x3_sample_unscaled = (max(X3) - min(X3)).*rand(N_points,1) + min(X3);
%x4_sample_unscaled = (max(X4) - min(X4)).*rand(N_points,1) + min(X4);

x1_sample_scaled = (x1_sample_unscaled-min(x1_sample_unscaled))/range(x1_sample_unscaled);
x2_sample_scaled = (x2_sample_unscaled-min(x2_sample_unscaled))/range(x2_sample_unscaled);
x3_sample_scaled = (x3_sample_unscaled-min(x3_sample_unscaled))/range(x3_sample_unscaled);
x4_sample_scaled = (x4_sample_unscaled-min(x4_sample_unscaled))/range(x4_sample_unscaled);

x_sample_scaled = [x1_sample_scaled x2_sample_scaled x3_sample_scaled x4_sample_scaled];
p_B_sample = modelfun(b_est, x_sample_scaled);

%%%%%%%%%%% first set of critical values %%%%%%%%%%%
k_L = find(p_B_sample<0.3);
p_B_sample_L = p_B_sample(k_L);
critical_T_L = x1_sample_unscaled(k_L);
critical_Fc_L = x2_sample_unscaled(k_L);
critical_Tc_L = x3_sample_unscaled(k_L);
critical_C_L = x4_sample_unscaled(k_L);

k_LL = find(p_B_sample>=0.3 & p_B_sample<0.45);
p_B_sample_LL = p_B_sample(k_LL);
critical_T_LL = x1_sample_unscaled(k_LL);
critical_Fc_LL = x2_sample_unscaled(k_LL);
critical_Tc_LL = x3_sample_unscaled(k_LL);
critical_C_LL = x4_sample_unscaled(k_LL);


k_ESD = find(p_B_sample>=0.45 & p_B_sample<0.85);
p_B_sample_ESD = p_B_sample(k_ESD);
critical_T_ESD = x1_sample_unscaled(k_ESD);
critical_Fc_ESD = x2_sample_unscaled(k_ESD);
critical_Tc_ESD = x3_sample_unscaled(k_ESD);
critical_C_ESD = x4_sample_unscaled(k_ESD);


%%%%%%%%%%%%%% Visualize histograms for each of the 4 levels
%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%% LOW LEVEL ALARM THRESHOLDS
%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
figure(6)
tiledlayout(2,2)
nexttile
histogram(critical_T_L)
xlabel('Temperature (K)')
ylabel('Occurences')
title('Histogram of Temperatures for which, p_B lies between 0.2 and 0.35')
subtitle('Critical Range: (501, 531]')

nexttile
histogram(critical_Fc_L)
xlabel('CW flow-rate (m^3/min)')
ylabel('Occurences')
title('Histogram of CW flow-rates for which, p_B lies between 0.2 and 0.35')
subtitle('Critical Range: (44.05, 44.64]')

nexttile
histogram(critical_Tc_L)
xlabel('CW temperature (K)')
ylabel('Occurences')
title('Histogram of CW temperatures for which, p_B lies between 0.2 and 0.35')
subtitle('Critical Range: (303.4, 303.86]')

nexttile
histogram(critical_C_L)
xlabel('Concentration (kmol/m^3)')
ylabel('Occurences')
title('Histogram of Concentrations for which, p_B lies between 0.2 and 0.35')
subtitle('Critical Range: [1.416, 1.554)')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERMEDIATE ALARM THRESHOLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%
figure(7)
tiledlayout(2,2)
nexttile
histogram(critical_T_LL)
xlabel('Temperature (K)')
ylabel('Occurences')
title('Histogram of Temperatures for which, p_B lies between 0.35 and 0.6')
subtitle('Critical Range: (472, 501]')

nexttile
histogram(critical_Fc_LL)
xlabel('CW flow-rate (m^3/min)')
ylabel('Occurences')
title('Histogram of CW flow-rates for which, p_B lies between 0.35 and 0.6')
subtitle('Critical Range: (43.47, 44.05]')

nexttile
histogram(critical_Tc_LL)
xlabel('CW temperature (K)')
ylabel('Occurences')
title('Histogram of CW temperatures for which, p_B lies between 0.35 and 0.6')
subtitle('Critical Range: (302.95, 303.4]')

nexttile
histogram(critical_C_LL)
xlabel('Concentration (kmol/m^3)')
ylabel('Occurences')
title('Histogram of Concentrations for which, p_B lies between 0.35 and 0.6')
subtitle('Critical Range: [1.554, 1.686)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIGH ALARM THRESHOLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%
figure(8)
tiledlayout(2,2)
nexttile
histogram(critical_T_LLL)
xlabel('Temperature (K)')
ylabel('Occurences')
title('Histogram of Temperatures for which, p_B lies between 0.6 and 0.75')
subtitle('Critical Range: (460, 472]')

nexttile
histogram(critical_Fc_LLL)
xlabel('CW flow-rate (m^3/min)')
ylabel('Occurences')
title('Histogram of CW flow-rates for which, p_B lies between 0.6 and 0.75')
subtitle('Critical Range: (43.22, 43.47]')

nexttile
histogram(critical_Tc_LLL)
xlabel('CW temperature (K)')
ylabel('Occurences')
title('Histogram of CW temperatures for which, p_B lies between 0.6 and 0.75')
subtitle('Critical Range: (302.76, 302.95]')


nexttile
histogram(critical_C_LLL)
xlabel('Concentration (kmol/m^3)')
ylabel('Occurences')
title('Histogram of Concentrations for which, p_B lies between 0.6 and 0.75')
subtitle('Critical Range: [1.686, 1.742]')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESD ALARM THRESHOLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%
figure(9)
tiledlayout(2,2)
nexttile
histogram(critical_T_ESD)
xlabel('Temperature (K)')
ylabel('Occurences')
title('Histogram of Temperatures for which, p_B is greater than 0.75')
subtitle('Critical Range: (458, 470]')

nexttile
histogram(critical_Fc_ESD)
xlabel('CW flow-rate (m^3/min)')
ylabel('Occurences')
title('Histogram of CW flow-rates for which, p_B is greater than 0.75')

nexttile
histogram(critical_Tc_ESD)
xlabel('CW temperature (K)')
ylabel('Occurences')
title('Histogram of CW temperatures for which, p_B is greater than 0.75')

nexttile
histogram(critical_C_ESD)
xlabel('Concentration (kmol/m^3)')
ylabel('Occurences')
title('Histogram of Concentrations for which, p_B is greater than 0.75')

%%%%% Creat histogram of Ts and Cs variables for p_B ~ 0.5 
k = find(p_B_sample>=0.45 & p_B_sample < 0.53);
figure(10)
tiledlayout(1,2)
nexttile
histogram(x1_sample_unscaled(k))
xlabel('Temperature (K)')
ylabel('Occurences')
title('Histogram of Temperatures for which, p_B ~ 0.5')

nexttile
histogram(x4_sample_unscaled(k))
xlabel('Concentration (kmol/m^3)')
ylabel('Occurences')
title('Histogram of Concentrations for which, p_B ~ 0.5')



%%%%%%%%%%%%% Regular data plots %%%%%%
figure(1)
tiledlayout(2,2)

nexttile
plot(X1, ydata, 'ko')
xlabel('Temperature (K)')
ylabel('Committer Probabilty, p_B')
title('p_B data as function of Temperature')

nexttile
plot(X2, ydata, 'ko')
xlabel('F_C (m^3/min)')
ylabel('Committer Probabilty, p_B')
title('p_B data as function of CW flow-rate')

nexttile
plot(X3, ydata, 'ko')
xlabel('T_C (K)')
ylabel('Committer Probabilty, p_B')
title('p_B data as function of CW temperature')

nexttile
plot(X4, ydata, 'ko')
xlabel('C_A (kmol/m^3)')
ylabel('Committer Probabilty, p_B')
title('p_B data as function of Concentration')

%%% Creat visual plots showing three different clusters and the alarm
%%% threshold demarcation line

figure(1)
plot(critical_T_L, p_B_sample_L, 'k*')
hold on
plot(critical_T_LL, p_B_sample_LL, 'r*')
hold on
plot(critical_T_ESD, p_B_sample_ESD, 'b*')
hold on
plot([509.5, 509.5], [0, 0.9], 'k--', 'LineWidth', 1)
hold on
plot([487.5, 487.5], [0, 0.9], 'r--', 'LineWidth', 1)
hold on
plot([453.1, 453.1], [0, 0.9], 'b--', 'LineWidth', 1)
xlabel('Reactor Temperature, T')
ylabel('Committer Probability, p_B')
title('Model predicted p_B as function of T')
legend('p_B < 0.3', '0.3 <= p_B < 0.45', '0.45 <= p_B < 0.85', 'Location','northeast')

figure(2)
plot(critical_Fc_L, p_B_sample_L, 'k*')
hold on
plot(critical_Fc_LL, p_B_sample_LL, 'r*')
hold on
plot(critical_Fc_ESD, p_B_sample_ESD, 'b*')
hold on
plot([44.21, 44.21], [0, 0.9], 'k--', 'LineWidth', 1)
hold on
plot([43.78, 43.78], [0, 0.9], 'r--', 'LineWidth', 1)
hold on
plot([43.1, 43.1], [0, 0.9], 'b--', 'LineWidth', 1)
xlabel('CW Flow-rate, F_c')
ylabel('Committer Probability, p_B')
title('Model predicted p_B as function of F_c')
legend('p_B < 0.3', '0.3 <= p_B < 0.45', '0.45 <= p_B < 0.85', 'Location','northeast')

figure(3)
plot(critical_Tc_L, p_B_sample_L, 'k*')
hold on
plot(critical_Tc_LL, p_B_sample_LL, 'r*')
hold on
plot(critical_Tc_ESD, p_B_sample_ESD, 'b*')
hold on
plot([303.52, 303.52], [0, 0.9], 'k--', 'LineWidth', 1)
hold on
plot([303.19, 303.19], [0, 0.9], 'r--', 'LineWidth', 1)
hold on
plot([302.66, 302.66], [0, 0.9], 'b--', 'LineWidth', 1)
xlabel('CW Temperature, T_c')
ylabel('Committer Probability, p_B')
title('Model predicted p_B as function of T_c')
legend('p_B < 0.3', '0.3 <= p_B < 0.45', '0.45 <= p_B < 0.85', 'Location','northeast')

