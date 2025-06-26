%% ������ǻ��ʧгƵ�ʣ�
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
% ǻ�� P

Pc = V_mc^2/(2*R_mc);
% ���� P
Pb = I0*U0;

beta_opt = Pb/Pc+1;   % ƥ�����ϵ��

Qext = Q_mc/beta_opt; % ��Ʒ������

R_L = R_mc/(1+beta_opt);
Q_L  = R_L/720;

%% HALF����ǻ ���ϵ���ݶ��� �� = 4.9e4

% ������ǻ��ʧг�ǣ���������ĸ��ؽ������㣬���ָ��ؽ�Ϊ -15��
% fai_s = HALF.fais_nat;  % ע������Ϊ��Ȼͬ����λʱ�������ؽ�Ϊ��ʱ��Vg����Iͬ��
fai_s = HALF.fais_mc_whc;

fai = fai_s/pi*180;
theta_L = -20 /180*pi;
Y = 2*I0/(V_mc/R_L);

tan_psi = Y*cos(fai_s)+tan(theta_L)*(1+Y*sin(fai_s)); % ��ʽ��ȷ
psi = atan(tan_psi)/pi*180
detuning_MC = tan_psi * HALF.f_rf /(2*Q_L)          % detuning frequency