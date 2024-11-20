close all;
clc;
clear;

%% The seventh task

N = 1e5;
c = 1;
Xi = -c.*rand(1, N);
Eta = c.*rand(1, N);
figure('Name', 'Совместная плотность вероятности ksi и eta');
plot(Xi, Eta, '.');
ylabel('\eta');
xlabel('\xi');
pax = gca;
pax.XLim = [-1.5 1.5];
pax.YLim = [-1.5 1.5];
grid on

hold on;
rectangle('Position', [-c, 0, c, c], 'EdgeColor', 'r');

expectedKsi = mean(Xi);
expectedEta = mean(Eta);
covMatrix = cov(Xi, Eta);
covKsiEta = covMatrix(1, 2);
fprintf('E[Xi] = %.5f, E[Eta] = %.5f, cov(Xi,Eta) = %.5f \n', expectedKsi, expectedEta, covKsiEta);

figure('Name', 'Гистограмма Xi');
% Nbins = int16(1+3.322*log10(N));
histogram(Xi,'Normalization','pdf');

figure('Name', 'Гистограмма Eta');
histogram(Eta,'Normalization','pdf');
