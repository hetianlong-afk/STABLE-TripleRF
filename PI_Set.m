% PI structure KP(s+KI/s)
%% Main cavity
% PIMC & PIHC
PIMC.m    = 80; % 取平均间隔 5120 , if 5120/2  32 tap FIR
% 64 tap FIR 
PIMC.tap  = 8;                            % for PI_Control_copy64
PIMC.dIQ  = PIMC.m/PIMC.tap;   % 一组IQ间隔     % for PI_Control_copy64
PIMC.d    = 800;  % 延迟间隔,必须是大于0的整数 800

PIMC.piIndex =0;
PIMC.Integral=0; % initial value PI 积分值
PIMC.RL      = HALF.R_mc; % 主腔负载阻抗

PIMC.Ts      = PIMC.m * HALF.Tb;
% PI.Ts = PI.dIQ * HALF.Tb;
% PI.KP = 1; % PI 控制器的比例系数 0.2 1 2 
% PI.KI = HALF.w_rf/(2*HALF.Q_mc)*PI.Ts*PI.KP/64;  % PI 控制器的积分系数 10 50 100

PIMC.KP = 1; % PI 控制器的比例系数 0.2 1 2 
PIMC.KI = 1e-1;  % PI 控制器的积分系数 10 50 100
% 测试发现KI要尽量小，以避免激励束流纵向振荡
PIMC.Ig0 = HALF.Ig_mc_0;                  % 发射机电流矢量
PIMC.DIg = 0;
PIMC.Ig  = PIMC.Ig0;                        % 初始 Ig=Ig0
PIMC.Vg0 = HALF.Vg_mc_0;
PIMC.DI  = 0;
PIMC.Ig_track=[];
%PIMC.Ig_FBid=[];
PIMC.Ig_FBid=zeros(1,PIMC.d/PIMC.m);
PIMC.RoverQ       = HALF.R_mc/HALF.Q_mc;  % 主腔的R/Q，用于PI反馈
PIMC.Vrf_mc_track = zeros(1,PIMC.m);  % 保持m个buckets看到的腔压矢量

% RF phase modulation 参数
% mf 调制度
PIMC.RFModulation_mf = 0.002;
% wm 调制频率
PIMC.RFModulation_wm = 0.0012*2/HALF.T0*2*pi; %0.0025  0.002505
% rf frequency
PIMC.RFModulation_wf = HALF.w_rf;
PIMC.RFModulation_Tb = HALF.Tb;
PIMC.Detun_time      = HALF.det_angle_mc/HALF.w_rf;

%% Harmonic cavity
PIHC.m    = 800; % 取平均间隔 5120 , if 5120/2  32 tap FIR
% 64 tap FIR 
PIHC.tap  = 8;                             % for PI_Control_copy64
PIHC.dIQ  = PIHC.m/PIHC.tap;   % 一组IQ间隔 % for PI_Control_copy64
PIHC.d    = 800;               % 延迟间隔,必须是大于0的整数 800

PIHC.piIndex =0;
PIHC.Integral=0; % initial value PI 积分值
PIHC.RL      = HALF.R_hc; % 主腔负载阻抗

PIHC.Ts      = PIHC.m * HALF.Tb;
PIHC.KP  = 1; % PI 控制器的比例系数 0.2 1 2 
PIHC.KI  = 4e-1;  % PI 控制器的积分系数 10 50 100
PIHC.Ig0 = HALF.Ig_hc_0;                  % 发射机电流矢量
PIHC.DIg = 0;
PIHC.Ig  = PIHC.Ig0;                        % 初始 Ig=Ig0
PIHC.Vg0 = HALF.Vg_hc_0;
PIHC.DI  = 0;
PIHC.Ig_track=[];
PIHC.Ig_FBid=zeros(1,PIHC.d/PIHC.m);
PIHC.RoverQ       = HALF.R_hc/HALF.Q_hc;  % 主腔的R/Q，用于PI反馈
PIHC.Vrf_hc_track = zeros(1,PIHC.m);  % 保持m个buckets看到的腔压矢量
PIHC.n_hc = HALF.n_hc;
%% Harmonic cavity2
PIHC2.m    = 800; % 取平均间隔 5120 , if 5120/2  32 tap FIR
% 64 tap FIR 
PIHC2.tap  = 8;                             % for PI_Control_copy64
PIHC2.dIQ  = PIHC2.m/PIHC2.tap;   % 一组IQ间隔 % for PI_Control_copy64
PIHC2.d    = 800;               % 延迟间隔,必须是大于0的整数 800

PIHC2.piIndex =0;
PIHC2.Integral=0; % initial value PI 积分值
PIHC2.RL      = HALF.R_hc2; % 主腔负载阻抗

PIHC2.Ts      = PIHC2.m * HALF.Tb;
PIHC2.KP  = 0.5; % PI 控制器的比例系数 0.2 1 2 
PIHC2.KI  = 3e-1;  % PI 控制器的积分系数 10 50 100
PIHC2.Ig0 = HALF.Ig_hc2_0;                  % 发射机电流矢量
PIHC2.DIg = 0;
PIHC2.Ig  = PIHC2.Ig0;                        % 初始 Ig=Ig0
PIHC2.Vg0 = HALF.Vg_hc2_0;
PIHC2.DI  = 0;
PIHC2.Ig_track=[];
PIHC2.Ig_FBid=zeros(1,PIHC2.d/PIHC2.m);
PIHC2.RoverQ       = HALF.R_hc2/HALF.Q_hc2;  % 主腔的R/Q，用于PI反馈
PIHC2.Vrf_hc_track = zeros(1,PIHC2.m);  % 保持m个buckets看到的腔压矢量
PIHC2.n_hc = n_hc2;







