%% Parámetros
clear all;

m11 = 0.6;
m12 = 0.03
m21 = 1.9;
m22 = 0.5;

p11 = -60;
p12 = -26;
p21 = 81;
p22 = -60;

S11 = m11 *exp(1i*((p11*pi)/180));
S12 = m12 *exp(1i*((p12*pi)/180));
S21 = m21 *exp(1i*((p21*pi)/180));
S22 = m22 *exp(1i*((p22*pi)/180));

%% Estabilidad
Delta = S11*S22 - S12*S21;
Delta_abs = abs(Delta)
K = (1-abs(S11)^2 - abs(S22).^2 + abs(Delta).^2)/(2*abs(S12*S21));

MUs = (1-abs(S22).^2) / ((S11-abs(Delta)*conj(S22)) + (S21*S12));
MUl = (1-abs(S11).^2) / ((S11-abs(Delta)*conj(S11)) + (S21*S12));

%% Error
U = (abs(S12 * S21 * S11 * S22)) / ((1-abs(S11).^2) * (1-abs(S22).^2));
u1 = 1 / (1+U).^2;
u2 = 1 / (1-U).^2;
U1 = 10*log10(u1);
U2 = 10*log10(u2);

%% Cálculo de Rs y Rl
%% En transistor unilateral
Rin = S11
Rout = S22
% Con máxima ganancia:
Rs = conj(S11)
Rl = conj (S22)

%% En transistor bilateral
B1 = 1 + abs(S11).^2 - abs(S22).^2 - abs(Delta).^2;
C1 = S11 - Delta * conj(S22);
Rs = (B1 + sqrt(B1^2 - 4*abs(C1).^2)) / 2*C1;
Rs_m = abs(Rs);
Rs_ph = rad2deg(angle(Rs));

B2 = 1 + abs(S22)^2 - abs(S11)^2 - abs(Delta)^2
C2 = S22 - Delta*conj(S11)
Rl = (B2 + sqrt(B2^2 - 4*abs(C2)^2)) / (2*C2)
Rl_m = abs(Rl);
Rl_ph = rad2deg(angle(Rl));

%% Circulo de ruido
Ropt_m = 0.62;
Ropt_p = 100;
Ropt = Ropt_m *exp(1i*((Ropt_p*pi)/180));

%Rs = ;
Rn = 20;
Z0 = 50;
F = 10.^(2/10);
Fmin = 10.^(1.6/10);

%N = (abs(Rs-Ropt).^2) / (1-abs(Rs).^2)
N = ((F-Fmin)/((4*Rn)/Z0))*abs(1+Ropt).^2;
%N = 10*log10(N)

Cf = Ropt_m/(N+1);
Cf_m = abs(Cf)
Cf_ph = rad2deg(angle(Cf))
Rf = (sqrt(N*(N+1-abs(Ropt).^2))) / (N+1);

%% Máxima ganancia
%% Caso Transistor unilateral
Gs_max = 1/(1-abs(Rs).^2);
Gl_max = 1/(1-abs(Rl).^2);
G0_max = abs(S21).^2;

% Ganancias en dB.
Gs_max_db = 10*log10(Gs_max)
G0_max_db =10*log10(G0_max)
Gl_max_db = 10*log10(Gl_max)

Gtu_max = Gs_max * G0_max * Gl_max;
Gtu_max_db = Gs_max_db + G0_max_db + Gl_max_db;

%% Caso Transistor bilateral
Rs_m = abs(Rs);
Rl_m = abs(Rl);
Rs_ph = rad2deg(angle(Rs));
Rl_ph = rad2deg(angle(Rl));

Rs = Rs_m *exp(1i*((Rs_ph*pi)/180));
Rl = Rl_m *exp(1i*((Rl_ph*pi)/180));

Gs_max = 1/(1-abs(Rs).^2);
G0_max = abs(S21).^2;
Gl_max = (1-abs(Rl).^2)/(abs(1-S22*Rl).^2);

% Ganancias en dB.
Gs_max_db = 10*log10(Gs_max)
G0_max_db =10*log10(G0_max)
Gl_max_db = 10*log10(Gl_max)

Gt_max = Gs_max * G0_max * Gl_max;
Gt_max_db = Gs_max_db + G0_max_db + Gl_max_db;

%% Circunferencias de ganancia máxima.
Gs = 10^(1/10);

gs = Gs / Gs_max

Cs = (gs*conj(S11))/(1-(1-gs)*abs(S11).^2)
Cs_m = abs(Cs)
Cs_ph = rad2deg(angle(Cs))
rs = (sqrt(1-gs)*(1-abs(S11).^2))/(1-(1-gs)*abs(S11).^2)