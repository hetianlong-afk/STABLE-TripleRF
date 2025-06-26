%% 计算主腔的失谐频率：
% Q_mc = 2.0e9;
% R_mc = Q_mc*90;
% I0     = 350e-3; 
% E0     = 2.2e9;  
% U0     = 600e3;  
% V_mc   = 1.4e6 ;  
clc;clear;
Q_mc = 2.0e9;
R_mc = Q_mc*720;
I0     = 200e-3; 
E0     = 5.5e9;  
U0     = 3.67e6;  
V_mc   = 11.3e6 ; 
% 腔耗 P

Pc = V_mc^2/(2*R_mc);
% 束耗 P
Pb = I0*U0;

beta_opt = Pb/Pc+1;   % 匹配耦合系数

Qext = Q_mc/beta_opt; % 外品质因素

R_L = R_mc/(1+beta_opt);
Q_L  = R_L/720;

%% HALF的主腔 耦合系数暂定在 β = 4.9e4

% 计算主腔的失谐角：根据所需的负载角来计算，保持负载角为 -15度
% fai_s = HALF.fais_nat;  % 注意能损；为自然同步相位时，当负载角为零时，Vg才与I同向。
fai_s = HALF.fais_mc_whc;

fai = fai_s/pi*180;
theta_L = -20 /180*pi;
Y = 2*I0/(V_mc/R_L);

tan_psi = Y*cos(fai_s)+tan(theta_L)*(1+Y*sin(fai_s)); % 此式正确
psi = atan(tan_psi)/pi*180
detuning_MC = tan_psi * HALF.f_rf /(2*Q_L)          % detuning frequency