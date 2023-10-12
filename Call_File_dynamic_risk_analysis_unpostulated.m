clc
clearvars
tau_init = 0.53;
N_runs = 200;
N_sim = 10;
output_var = zeros(N_runs,5);
output_L_failure_probability = zeros(N_runs,1);
output_ESD_failure_probability = zeros(N_runs,1);

tic
parpool('local',4);

parfor i = 1:3

    [N_L, N_ESD, N_L_to_ESD, N_basin_B, time_L_ESD, time_ESD_basin_B, Avg_N_B ] = dynamic_risk_analysis_unpostulated(N_sim, tau_init);
    output_var(i, :) = [N_L, N_ESD, N_L_to_ESD, N_basin_B, Avg_N_B ];
    output_L_failure_probability(i) = N_L_to_ESD/N_L;
    output_ESD_failure_probability(i) = N_basin_B/N_ESD;
end
toc
delete(gcp('nocreate'));

m = find(output_L_failure_probability>0.16 & output_L_failure_probability<=0.26);
pd = fitdist(output_L_failure_probability(m), 'Beta');
y = pdf(pd, x);



x_data = output_L_failure_probability(m);
figure(1)
h = histogram(output_L_failure_probability(m))
xlabel('L Failure Probabiltity, p_{failure, L}')
ylabel('Number of Occurrences')
title('Histogram of 150 L Failure Probabilities over 15000 simulations')
figure(3)

plot(x, y)



figure(2)
h1 = histfit(output_L_failure_probability(m), 9, "beta")
%hold on
%plot(x,y)
legend('Data', 'Beta PDF fit ')
xlabel('L Failure Probabiltity, p_{failure, L}')
ylabel('Number of Occurrences')


figure(3)
histogram(x_data,9)
hold on
plot(x, y)
xlim([0.13, 0.3])

x = 0:100;
y = binopdf(x, 100, 0.002);
figure(3)
plot(x./100,y)