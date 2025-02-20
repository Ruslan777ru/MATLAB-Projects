
clc;
close all;
clear;

%% Расчет смесителя

fprintf('\n ____________Расчет смесителя______________ \n')

Gs = 0.02;
GH = 0.02;

gamma = 27;
Ug = 406e-3;
I0 = 8e-9;

arg = gamma * Ug; % вычисляем аргумент для функции Бесселя
J0 = besseli(0, arg);
J1 = besseli(1, arg);

G0 = gamma * I0 * besseli(0, arg);
G1 = gamma * I0 * besseli(1, arg);

% G0 = gamma * I0 * J0;
% G1 = gamma * I0 * J1;


disp(['G0 = ', num2str(G0)])
disp(['G1 = ', num2str(G1)])

K = abs((4 * G1 * Gs^2 * GH) / ((G0 * (GH + Gs) + 2*GH*Gs)^2 ...
    - 0.5*G1^2*(GH + Gs)^2));

disp(['K = ', num2str(K)])

Kpr = real(20 * log10(K));
disp(['Kpr = ', num2str(Kpr)])


U_vh = 34.895e-3;
U_vih = 1.987e-3;

K1 = U_vih/U_vh

K_real1 = 20*log10(K1)

K_sravnenie1 = K / K1

K_sravn_dB1 = 20 * log10(K_sravnenie1)
