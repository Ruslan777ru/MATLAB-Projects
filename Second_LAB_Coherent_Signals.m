clc;
clear;
close all;
%%

fs = 1000; %Частота дискретизации
fc = 100; %Частота несущей
T = 1; %Длительность импульса
N = T * fs; %Количество отсчетов
f_offset = 0; % Нулевой сдвиг по частоте
t = linspace(0, T, N); %Создание когерентной пачки радиоимпульсов
rect_pulse = rectpuls(t, T); %Прямоугольный импульс
t_c = 0:(1/fs):(10-1/fs);
signal_coherent = cos(2*pi*fc*t);
signal_coherent_with_rect = signal_coherent .* rect_pulse;
signal_coherent_with_offset = signal_coherent_with_rect .* exp(1i*2*pi*f_offset*t);

signal_incoherent = randn(1, N); %Создание некогерентной пачки радиоимпульсов
signal_incoherent_with_rect = signal_incoherent .* rect_pulse;


LEN_coherent = length(signal_coherent); %Создание АЧХ для когернетного сигнала
FFT_coherent = fft(signal_coherent);
FREQ_coherent = (0:LEN_coherent-1)*(fs/LEN_coherent);
POWER_coherent = abs(FFT_coherent).^2/LEN_coherent;

LEN_incoherent = length(signal_incoherent); %Создание АЧХ для некогернетного сигнала
FFT_incoherent = fft(signal_incoherent);
FREQ_incoherent = (0:LEN_incoherent-1)*(fs/LEN_incoherent);
POWER_incoherent = abs(FFT_incoherent).^2/LEN_incoherent;

figure; 
subplot(2,1,1);
plot(FREQ_coherent, POWER_coherent);
title('АЧХ когерентного сигнала');
xlabel('Частота');
ylabel('Энергия');
grid on

subplot(2,1,2);
plot(FREQ_incoherent, POWER_incoherent);
title('АЧХ некогерентного сигнала');
xlabel('Частота');
ylabel('Энергия');
grid on


figure;
subplot(2, 1, 1);
plot(t, signal_coherent_with_rect);
title('Когерентная пачка радиоимпульсов');
xlabel('Время');
ylabel('Амплитуда');
grid on


subplot(2, 1, 2);
plot(t, signal_incoherent_with_rect);
title('Некогерентная пачка радиоимпульсов');
xlabel('Время');
ylabel('Амплитуда');
grid on

autocorr_coherent = xcorr(signal_coherent_with_rect, 'biased'); %Расчет и построение тела неопределенности и его сечений
autocorr_incoherent = xcorr(signal_incoherent_with_rect, 'biased');

figure; %Сечение при нулевом сдвиге по времени
subplot(2, 1, 1);
plot(t, abs(autocorr_coherent(N:2*N-1)));
title('Сечение при нулевом сдвиге по времени (когерентная)');
xlabel('Время');
ylabel('Амплитуда');
grid on

subplot(2, 1, 2);
plot(t, abs(autocorr_incoherent(N:2*N-1)));
title('Сечение при нулевом сдвиге по времени (некогерентная)');
xlabel('Время');
ylabel('Амплитуда');
grid on

%Сдвиг при нулевом сдвиге по частоте

FFT_Shift_coherent = fftshift(FFT_coherent);
FREQ_Shift_coherent = (-LEN_coherent/2:LEN_coherent/2-1)*(fs/LEN_coherent);
POWER_Shift_coherent = abs(FFT_Shift_coherent).^2/LEN_coherent;

FFT_Shift_incoherent = fftshift(FFT_incoherent);
FREQ_Shift_incoherent = (-LEN_incoherent/2:LEN_incoherent/2-1)*(fs/LEN_incoherent);
POWER_Shift_incoherent = abs(FFT_Shift_incoherent).^2/LEN_incoherent;

figure;
subplot(2,1,1);
plot(FREQ_Shift_coherent, POWER_Shift_coherent);
title('Сечение при нулевом сдвиге по частоте (когерентная)');
xlabel('Частота');
ylabel('Энергия');
grid on;

subplot(2,1,2);
plot(FREQ_Shift_incoherent, POWER_Shift_incoherent);
title('Сечение при нулевом сдвиге по частоте (некогерентная)');
xlabel('Частота');
ylabel('Энергия');
grid on;

%Реализация обнаружителя на основе коррелятора
threshold = 0.5; %Порог обнаружения

correlation_coherent = xcorr(signal_coherent_with_rect, signal_coherent_with_rect, 'biased'); %Обнаружитель на основе коррелятора с опорным сигналом
detection_coherent = max(correlation_coherent) > threshold;

correlation_incoherent = xcorr(signal_incoherent_with_rect, signal_incoherent_with_rect, 'biased'); %Обнаружитель на основе коррелятора с принимаемым сигналом
detection_incoherent = max(correlation_incoherent) > threshold;

disp(['Обнаружение на основе коррелятора с опорным сигналом: ', num2str(detection_coherent)]);
disp(['Обнаружение на основе коррелятора с принимаемым сигналом: ', num2str(detection_incoherent)]);

mean_correlation_coherent = mean(correlation_coherent);
std_correlation_coherent = std(correlation_coherent);

mean_correlation_Necoherent = mean(correlation_incoherent);
std_correlation_Necoherent = std(correlation_incoherent);

Pfa_values = [1e-1, 1e-4, 1e-6];
thresholds = zeros(size(Pfa_values));
Pd_values_coherent = zeros(size(Pfa_values));
Pd_values_Necoherent = zeros(size(Pfa_values));

for i = 1:length(Pfa_values)
    thresholds(i) = mean_correlation_coherent + sqrt(2) * std_correlation_coherent * qfuncinv(Pfa_values(i));
    Pd_values_coherent(i) = qfunc((thresholds(i) - mean_correlation_coherent) / std_correlation_coherent);
end

disp('Кривая оценки вероятности ложной тревоги для когерентной пачки радиоимпульсов:');
disp([Pfa_values; thresholds]);


for i = 1:length(Pfa_values)
    thresholds(i) = mean_correlation_Necoherent + sqrt(2) * std_correlation_Necoherent * qfuncinv(Pfa_values(i));
    Pd_values_Necoherent(i) = qfunc((thresholds(i) - mean_correlation_Necoherent) / std_correlation_Necoherent);
end

disp('Кривая оценки вероятности правильного обнаружения для когерентной пачки радиоимпульсов:');
disp([Pfa_values; Pd_values_coherent]);
disp([Pfa_values; Pd_values_Necoherent]);

