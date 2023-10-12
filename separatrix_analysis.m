clc

tau = 0.53;
tBounds = [0 50];
h = 0.01;
ts = tBounds(1):h:tBounds(2);
Fc_0 = 50;
Tc_0 = 300;
Kc = 0.02;


N_sim = 150000;

Cs_init = (1.99 - 0).*rand(N_sim,1) + 0;
Ts_init = (1000 - 300).*rand(N_sim,1) + 300;
Tc_init = (310 - 301).*rand(N_sim,1) + 301;
Fc_init = (52 - 39).*rand(N_sim,1) + 39;

Ts_ss = zeros(N_sim, 1);
Cs_ss = zeros(N_sim, 1);
Fc_ss = zeros(N_sim, 1);
Tc_ss = zeros(N_sim, 1);

[Cs, Ts, Fc, Tc] = CSTR_P_control(tau, tBounds, [1.1, 750, Fc_0, Tc_0], Fc_0, Tc_0, Kc);

tic
parpool('local', 10);
parfor k = 1:N_sim

    [Cs, Ts, Fc, Tc] = CSTR_P_control(tau, tBounds, [Cs_init(k), Ts_init(k), Fc_init(k), Tc_init(k)], Fc_0, Tc_0, Kc);
    Ts_ss(k) = Ts(end);
    Cs_ss(k) = Cs(end);
    Fc_ss(k) = Fc(end);
    Tc_ss(k) = Tc(end);

end
toc
delete(gcp('nocreate'));

k_A = find(Ts_ss >800);
k_B = find(Ts_ss<= 400);

responses = zeros(N_sim, 1);
responses(k_A) = 1;
responses(k_B) = 0;

%% randomly assign training and testing data %%
p_train = 0.8;
p_test = 1 - p_train;
k_train = transpose(randperm(N_sim, p_train*N_sim));
k_total = transpose(1:N_sim);
k_total(k_train) = [];
k_test = k_total;

%%%%%

data_SVM = zeros(N_sim, 5);
data_SVM (:, 1:4) = [Cs_init Ts_init Fc_init Tc_init];
data_SVM(:, 5) = responses;

data_train = data_SVM(k_train,:);
data_test = data_SVM(k_test, :);

tic
SVM_mdl = fitcsvm(data_train(:, 1:2), data_train(:,5));
toc

pred = predict(SVM_mdl, data_test(:, 1:2));



k= find(pred==data_test(:,5));
SVM_accuracy = 100*length(k)/(p_test*N_sim)

%% False positives and false negatives 
k_false = find(pred~=data_test(:,5));
pred_false = pred(k_false);
data_false = data_test(k_false,:);


k_false_positive = find(pred_false == 1);
k_false_negative = find(data_false(:,5) == 1);

data_false_positive = data_false(k_false_positive, :);
data_false_negative = data_false(k_false_negative, :);






Conc_sep = data_SVM(:,1);
Temp_sep = data_SVM(:,2);
Fc_init = data_SVM(:, 3);
Tc_init = data_SVM(:, 4);


%k = find(p_B_sample>=0.45 & p_B_sample < 0.53);
k1 = find(Conc>1.615 & Conc<1.656 & Temp>478 & Temp<488);
%k1 = find(Conc >= 1.636 & Conc <= 1.647 & Temp>= 480.6 & Temp <= 482.8);




figure(3)
gscatter(Conc_sep, Temp_sep, responses, 'br')
hold on
%plot(Conc(k1),Temp(k1), 'k*')
plot(Conc, Temp, 'k*')
hold on
legend('Basin B steady-states', 'Basin A steady-states', 'Model Values for which, p_B ~ 0.5', 'location','northeast')
xlabel('Initial values of Concentration (kmol/m^3)')
ylabel('Initial Values of Temperature (K)')

figure(4)
gscatter(Temp, Fc_init, responses, 'br')
hold on
plot(Temp(k1),Fc_init(k1), 'k*')
hold on
legend('Basin B steady-states', 'Basin A steady-states', 'Model Values for which, p_B ~ 0.5', 'location','northeast')
xlabel('Initial values of Concentration (kmol/m^3)')
ylabel('Initial Values of Temperature (K)')

Conc_pB_05 = Conc(k1);
Temp_pB_05 = Temp(k1);

%%%% Scatter plot of false positives and false negatives %%%
figure(4)
scatter(data_false_positive(:,1),data_false_positive(:,2), 'ro')
xlabel('Concentration (kmol/m^3)')
ylabel('Temperature (K)')
title('False Positive Variable Values')

figure(5)
scatter(data_false_negative(:,1),data_false_negative(:,2), 'ko')
xlabel('Concentration (kmol/m^3)')
ylabel('Temperature (K)')
title('False Negativee Variable Values')
