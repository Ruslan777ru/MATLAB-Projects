
clc;
clear;
close all;

%% Необходимые параметры транзистора 2Т903А

% Параметры идеализированных статических характеристик
f = 8e6;
f_min = 1.6e6;
f_max = 8e6;

% Сопротивления
r_nas = 1;
r_b = 2;
r_e = 0;
R_ue = 0.1;
h_21 = 15;

% Высокочастотные параметры
f_x     = 125e6;

C_k     = 120e-12;
C_ka    = 30e-12;
C_e     = 400e-12;
tau_k   = 500e-12;

L_e     = 5e-9;
L_b     = 20e-9;
L_k     = 5e-9;

% Допустимые параметры
E_kb    = 0;
E_ke    = 60;
E_k_dop = 0;
E_b_dop = 4;

I_ko_dop = 3;
I_bo_dop = 0;
I_k_max  = 10;

delta_F  = 0;

% Тепловые параметры
temp_dop = 150;
R_pk     = 3.33;

% Экспериментальные параметры
f_strih = 50e6;
P_k_strih = 450;
K_p = 15;
KPD_tr = 65;
E_k_strih = 50;

% _________ Схема включения: ОЭ
% _________ Режим работы: Класс В
%__________ Линейный <-30 дБ

% Расчет необходимой мощности
P_fider = 5;
KPD     = 0.9;
KBV     = 0.85;
KSV     = 1.18;
R_fider = 75;

% Угол отсечки (тк класс Б то можно взять 90)
teta = 90;

%____Коэффициенты Берга_____%
%
% t_korpus = 20; 
% teta_v = 90;
% teta_v = teta_v *pi/180;
% 
% Lambda0v = (sin(teta_v)-teta_v*cos(teta_v))/(pi*(1-cos(teta_v)))
% Lambda1v = (teta_v-sin(teta_v)*cos(teta_v))/(pi*(1-cos(teta_v)))
% 
%___________________________%

% Коэффициенты Берга
alpha_1 = 0.5;
gamma_1 = alpha_1;

alpha_0 = 0.32;
gamma_0 = alpha_0;

% Напряжение питания
Ek = 25;

%% Расчет
fprintf('\n ____________РАСЧЕТ ФИДЕРА______________ \n')

% Необходимая выходная мощность
P_out = P_fider / (KPD * KBV);
disp(['Мощность на фидере P_out: ' num2str(P_out)])

% Мощность на двухтактном усилителе на одном из двух транзисторов
P_one = P_out / 2;
disp(['Мощность на двухтактном усилителе на одном' ...
    ' из двух транзисторов P_one: ' num2str(P_one)])


%% Расчет коллекторной цепи
fprintf('\n ____________РАСЧЕТ КОЛЛЕКТОРНОЙ ЦЕПИ______________ \n')

% Амплитуда напряжений первой гармоники на коллекторе:
U_k1_gr = Ek* (0.5 + 0.5 * ...
    sqrt (1 - ( (8 * r_nas * P_one) / (alpha_1 * Ek.^2) ) ));

disp(['Амплитуда напряжений первой гармоники на коллекторе U_k1_gr: ' ...
    num2str(U_k1_gr)])

% Максимальное напряжение на коллекторе
E_k_max = Ek + 1.3 * U_k1_gr;
disp(['Максимальное напряжение на коллекторе: ' num2str(E_k_max)])

if E_k_max < E_ke
    disp(['E_k_max мал по сравнению с E_ke. ПРАВИЛЬНО!'])
else 
     disp(['E_k_max большой по сравнению с E_ke. ОШИБКА!!!!! '])
end

% Амплитуда первой гармоники тока
fprintf('\n')
I_k1 = 2 * P_one / U_k1_gr;
disp(['Амплитуда первой гармоники тока I_k1: ' num2str(I_k1)])


% Постоянная составляющая коллекторного тока
I_k0 = I_k1 * alpha_0 / alpha_1;
disp(['Постоянная составляющая коллекторного тока I_k0: ' num2str(I_k0)])

if I_k0 < I_ko_dop
     disp(['I_k0 мал по сравнению с I_ko_dop. ПРАВИЛЬНО!'])
else 
     disp(['I_k0 большой по сравнению с I_ko_dop. ОШИБКА!!!!! '])
end


% Максимальный коллекторный ток 
fprintf('\n')
I_k_max_col = I_k0 / alpha_0;
disp(['Максимальный коллекторный ток I_k_max_col : ' num2str(I_k_max_col )])

if I_k_max_col < I_k_max
    disp(['I_k_max_col мал по сравнению с I_k_max. ПРАВИЛЬНО!'])
else 
    disp(['I_k_max_col большой по сравнению с I_k_max. ОШИБКА!!!!! '])
end

% Максимальная мощность, потребляемая от источника питания
fprintf('\n')
P_0_max = Ek * I_k0;
disp(['Максимальная мощность, потребляемая от источника питания' ...
    ' P_0_max : ' num2str(P_0_max )])

% КПД коллекторной цепи
KPD_col = P_one / P_0_max;
disp(['КПД коллекторной цепи KPD_col : ' num2str(KPD_col)])

% Максимальная рассеиваемая мощность на коллекторе транзистора
P_k_max = P_0_max - P_one * sqrt(KBV);
disp(['Максимальная рассеиваемая мощность на коллекторе транзистора P_k_max : ' ...
    num2str(P_k_max)])

% Номинальное сопротивление коллекторной нагрузки в одном из плеч
R_ek_nom = U_k1_gr.^2 / (2 * P_one);
disp(['Номинальное сопротивление коллекторной нагрузки в одном из плеч: ' ...
    num2str(R_ek_nom)])

R_n = 2 * R_ek_nom;
disp(['Cопротивление коллекторной нагрузки R_n: ' ...
    num2str(R_n)])

%% Расчет входной цепи
fprintf('\n ____________РАСЧЕТ ВХОДНОЙ ЦЕПИ______________ \n')

% Сопротивление R2 
R2 = h_21 / (2 * pi* f_x * C_e) * (1 - (h_21 / (2 * pi * f_x * R_ue) ) );
disp(['Сопротивление R2: ' num2str(R2)])

% Сопротивление R1
R1 = h_21 / (2 * pi* f_x * C_k);
disp(['Сопротивление R1: ' num2str(R1)])

% Расчет параметра Каппа
Cappa = 1 + gamma_1 * 2 * pi * f_x * C_k * R_ek_nom;
disp(['Расчет параметра Каппа: ' num2str(Cappa)])

% Амплитуда тока базы на частоте 8 МГц
I_b = Cappa * (sqrt(1 + (h_21 * f/f_x).^2) / h_21 * gamma_1) * I_k1; 
disp(['Амплитуда тока базы I_b: ' num2str(I_b)])

% Расчет параметров L_вхОЭ, r_вхОЭ, С_вхОЭ, r_э = 0
fprintf('\n')
fprintf('Расчет параметров L_вхОЭ, r_вхОЭ, С_вхОЭ, r_э = 0\n')
L_vhOE = L_b + L_e/Cappa;
r_vhOE = ( (1 + gamma_1*2*pi*f_x*C_ka*R_ek_nom) * r_b + r_e + gamma_1 * 2 * pi * f_x * L_e) / Cappa;
R_vhOE = (r_b + (1 + gamma_1 * h_21) * r_e) / Cappa - r_vhOE + R2 * (1 - gamma_1);
C_vhOE = h_21 / (2 * pi * f_x * R_vhOE);

disp(['Расчет параметра L_vhOE: ' num2str(L_vhOE)])
disp(['Расчет параметра r_vhOE: ' num2str(r_vhOE)])
disp(['Расчет параметра R_vhOE: ' num2str(R_vhOE)])
disp(['Расчет параметра C_vhOE: ' num2str(C_vhOE)])

% Реактивная и резистивная составляющие сопротивления
fprintf('\n')
fprintf('Реактивная и резистивная составляющие сопротивления\n')
r_vh     = r_vhOE + R_vhOE / (1 + (h_21 * f/f_x).^2);
Cappa_vh = 2*pi*f*L_vhOE + (R_vhOE*(h_21 * f/f_x) / (1 + (h_21*f/f_x)).^2);

disp(['Расчет параметра r_vh: ' num2str(r_vh)])
disp(['Расчет параметра Cappa_vh: ' num2str(Cappa_vh)])

% Входная мощность
P_vh = 0.5 * I_b.^2 * r_vh;
disp(['Входная мощность P_vh: ' num2str(P_vh)])

% Коэффициент усиления по мощности
Kp = P_one / P_vh; 
disp(['Коэффициент усиления по мощности Kp: ' num2str(Kp)])

% Постоянное составляющие базового и эмиттерного токов
fprintf('\n')
fprintf('Постоянное составляющие базового и эмиттерного токов\n')
I_b0 = I_k0 / h_21;
I_e0 = I_b0 + I_k0;

disp(['Расчет параметра I_b0: ' num2str(I_b0)])
disp(['Расчет параметра I_e0: ' num2str(I_e0)])

% Максимальная мощность, рассеиваемая в транзисторе
fprintf('\n')
P_ras = P_vh + P_k_max;
disp(['Максимальная мощность, рассеиваемая в ' ...
    'транзисторе P_ras: ' num2str(P_ras)])

% Блокировочный ддроссель на источник коллекторного питания
L_dr = 10 * R_n / (2 * pi * f_min);
disp(['L_dr >= ' num2str(L_dr) '      L_dr = 0.001'])
L_dr = 0.001;

% Конденсатор C_bl на источнике коллекторного питания
C_bl = 10 / (2 * pi * f_min * R_n);
disp(['C_bl >= ' num2str(C_bl) '     C_bl = 8e-9'])
C_bl = 8e-9;

% Конденсатор C_p1, блокирующий протекание постоянной составляющей
% предыдущих каскадов
C_p1 = 10 / (2 * pi * f_min * R_vhOE);
disp(['C_p1 >= ' num2str(C_p1) '     C_p1 = 6e-9'])
C_p1 = 6e-9;

%% Проектирование выходной колебательной системы
fprintf('\n_______ПРОЕКТИРОВАНИЕ ВЫХОДНОЙ КОЛЕБАТЕЛЬНОЙ СИСТЕМЫ________\n')

% Коэффициент перекрытия 
Kf = f_max / f_min;
disp(['Коэффициент перекрытия Kf: ' num2str(Kf)])
fprintf('\n')

R1_vih = R_n;     % выходное сопротивление
R2_vih = R_fider; % сопротивление нагрузки

% Для транзистора 2Т903А 
% Коэффициент перекрытия 1.5. Полученные промежутки:
% 1.6 - 2.4, 2.4 - 3.5, 3.5 - 5.3, 5.3 - 8

% Промежуточные сопротивления на частоту 8 МГц
R_p1 = 10;
R_p2 = 35;
R_p3 = 45;

disp(['Промежуточное сопротивление Rп1: ' num2str(R_p1)])
disp(['Промежуточное сопротивление Rп2: ' num2str(R_p2)])
disp(['Промежуточное сопротивление Rп3: ' num2str(R_p3)])

% Промежуточные сопротивления трансформации внутри Т-цепочек

%____Выбраны рандомно____%
R_01 = 200;
R_02 = 400;
R_03 = 600;
R_04 = 800;
%________________________%

% Первый промежуток сопротивлений
X_L1 = R1_vih * sqrt(R_01/R1_vih - 1);
X_L2 = R_p1   * sqrt(R_01/R_p1 - 1);
X_C1 = R_01 / ( sqrt(R_01/R1_vih - 1) + sqrt(R_01/R_p1 - 1));

% Второй промежуток сопротивлений
X_L3 = R_p1   * sqrt(R_02/R_p1 - 1);
X_L4 = R_p2   * sqrt(R_02/R_p2 - 1);
X_C2 = R_02 / ( sqrt(R_02/R_p1 - 1) + sqrt(R_02/R_p2 - 1));

% Третий промежуток сопротивлений
X_L5 = R_p2   * sqrt(R_03/R_p2 - 1);
X_L6 = R_p3   * sqrt(R_03/R_p3 - 1);
X_C3 = R_03 / ( sqrt(R_03/R_p2 - 1) + sqrt(R_03/R_p3 - 1));

% Четвертый промежуток сопротивлений
X_L7 = R_p3   * sqrt(R_04/R_p3 - 1);
X_L8 = R2_vih   * sqrt(R_04/R2_vih - 1);
X_C4 = R_04 / ( sqrt(R_04/R_p3 - 1) + sqrt(R_04/R2_vih - 1));

fprintf('\nПервый промежуток сопротивлений\n')
disp(['X_L1: ' num2str(X_L1)])
disp(['X_L2: ' num2str(X_L2)])
disp(['X_C1: ' num2str(X_C1)])

fprintf('\nВторой промежуток сопротивлений\n')
disp(['X_L3: ' num2str(X_L3)])
disp(['X_L4: ' num2str(X_L4)])
disp(['X_C2: ' num2str(X_C2)])

fprintf('\nТретий промежуток сопротивлений\n')
disp(['X_L5: ' num2str(X_L5)])
disp(['X_L6: ' num2str(X_L6)])
disp(['X_C3: ' num2str(X_C3)])

fprintf('\nЧетвертый промежуток сопротивлений\n')
disp(['X_L7: ' num2str(X_L7)])
disp(['X_L8: ' num2str(X_L8)])
disp(['X_C4: ' num2str(X_C4)])


% Расчет компонентной базы
L_1 = X_L1 / (2*pi*f);
L_2 = X_L2 / (2*pi*f);
C_1 = 1 / (2*pi*f*X_C1);

L_3 = X_L3 / (2*pi*f);
L_4 = X_L4 / (2*pi*f);
C_2 = 1 / (2*pi*f*X_C2);

L_5 = X_L5 / (2*pi*f);
L_6 = X_L6 / (2*pi*f);
C_3 = 1 / (2*pi*f*X_C3);

L_7 = X_L7 / (2*pi*f);
L_8 = X_L8 / (2*pi*f);
C_4 = 1 / (2*pi*f*X_C4);

L_mid1 = L_2 + L_3;
L_mid2 = L_4 + L_5;
L_mid3 = L_6 + L_7;

fprintf('\nПервый промежуток\n')
disp(['L_1: ' num2str(L_1)])
disp(['L_2: ' num2str(L_2)])
disp(['C_1: ' num2str(C_1)])

fprintf('\nВторой промежуток\n')
disp(['L_3: ' num2str(L_3)])
disp(['L_4: ' num2str(L_4)])
disp(['C_2: ' num2str(C_2)])

fprintf('\nТретий промежуток\n')
disp(['L_5: ' num2str(L_5)])
disp(['L_6: ' num2str(L_6)])
disp(['C_3: ' num2str(C_3)])

fprintf('\nЧетвертый промежуток\n')
disp(['L_7: ' num2str(L_7)])
disp(['L_8: ' num2str(L_8)])
disp(['C_4: ' num2str(C_4)])

fprintf('\nЗначения по середине согласованной цепи\n')
disp(['L_mid1: ' num2str(L_mid1)])
disp(['L_mid2: ' num2str(L_mid2)])
disp(['L_mid3: ' num2str(L_mid3)])


% Заголовки для согласованной цепи в Microcap 12
fprintf('\n_______ЗАГОЛОВКИ ДЛЯ СОГЛАСОВАННОЙ ЦЕПИ В MICROCAP 12________\n')
disp(['.define R1_vih ' num2str(R1_vih)])
disp(['.define R2_vih ' num2str(R2_vih)])

disp(['.define L_1 ' num2str(L_1)])
disp(['.define L_2 ' num2str(L_2)])
disp(['.define C_1 ' num2str(C_1)])

disp(['.define L_3 ' num2str(L_3)])
disp(['.define L_4 ' num2str(L_4)])
disp(['.define C_2 ' num2str(C_2)])

disp(['.define L_5 ' num2str(L_5)])
disp(['.define L_6 ' num2str(L_6)])
disp(['.define C_3 ' num2str(C_3)])

disp(['.define L_7 ' num2str(L_7)])
disp(['.define L_8 ' num2str(L_8)])
disp(['.define C_4 ' num2str(C_4)])

disp(['.define L_mid1 ' num2str(L_mid1)])
disp(['.define L_mid2 ' num2str(L_mid2)])
disp(['.define L_mid3 ' num2str(L_mid3)])