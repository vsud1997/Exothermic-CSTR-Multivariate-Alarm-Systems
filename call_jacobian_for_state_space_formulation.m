clc
clearvars
tau = 0.5;
tBounds = [0 50];
h = 0.01;
ts = tBounds(1):h:tBounds(2);
Fc_0 = 50;
Tc_0 = 300;
Kc = 0.02;
k0 = 17.038; % Rate constant
E = 1.50e+04; % Actication energy
R = 8.314; % Gas Constant
alpha = 0.075; 
beta = 9370.9;
V_j = 10;

T_sp = 800;

[Cs, Ts, Fc, Tc] = CSTR_P_control(tau, tBounds, [1.1, 750, Fc_0, Tc_0], Fc_0, Tc_0, Kc);

C_ss = Cs(end);
T_ss = Ts(end);
Fc_ss = Fc(end);
Tc_ss = Tc(end);

[A,B] = jacobian_CSTR_P_controlled(k0, tau, Kc, Fc_0, Tc_0,R, alpha, beta, V_j, E, tBounds, C_ss, T_ss, Fc_ss, Tc_ss);

N_sim = 10000;

X = zeros(3,N_sim);
U = zeros(1, N_sim - 1);

x_init = [1.1; 750; 300];

X(:, 1) = x_init;
U(1) = Fc_0;

for i = 2:N_sim

    
end

figure(1)
plot(1:N_sim, X(2,:), 'ko')

