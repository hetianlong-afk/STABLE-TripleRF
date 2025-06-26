function [HALF]=machine(C,I0,U0,E0,tau_s,tau_z,sigma_t0,sigma_e0,alpha_c,h,V_mc,V_hc,V_hc2,fais_mc,fais_hc,fais_hc2,n_hc,n_hc2,R_hc,Q_hc,R_hc2,Q_hc2,fillrate,fre_shift_hc,fre_shift_hc2,Q_mc,R_mc,fre_shift_mc,F1,F2,F0)
% ��������������������HALF�ṹ����
cspeed = 299792458;
HALF.C       = C;
HALF.R       = C/(2*pi);
HALF.I0      = I0;
HALF.h       = h;
HALF.U0      = U0;
HALF.E0      = E0;
HALF.sigma_t0= sigma_t0;
HALF.sigma_e0= sigma_e0;
HALF.alpha_c = alpha_c;
HALF.tau_s   = tau_s;
HALF.tau_z   = tau_z;
HALF.V_mc    = V_mc;
HALF.V_hc    = V_hc;
HALF.V_hc2   = V_hc2;
HALF.fais_mc = fais_mc;
HALF.fais_hc = fais_hc;
HALF.fais_hc2= fais_hc2;
HALF.n_hc    = n_hc;
HALF.n_hc2   = n_hc2;

HALF.R_hc  = R_hc;  
HALF.Q_hc   = Q_hc;
HALF.R_hc2   = R_hc2;  
HALF.Q_hc2   = Q_hc2;
HALF.fillrate  = fillrate;
HALF.fre_shift_hc = fre_shift_hc;
HALF.fre_shift_hc2 = fre_shift_hc2;
HALF.R_mc    = R_mc;  
HALF.Q_mc    = Q_mc;
HALF.fre_shift_mc = fre_shift_mc;

HALF.gamma = E0/0.51099906/1e6;       % ���Ӿ�ֹ����  0.51099906 MeV
HALF.beta  = sqrt(1-1/(HALF.gamma^2));

HALF.T0 = HALF.C/cspeed/HALF.beta;    % ��������
HALF.Tb = HALF.T0/HALF.h;             % ���� norminal bucket position ����
HALF.f_rf = HALF.h/HALF.T0;           % ��ƵƵ��
HALF.w_rf = HALF.f_rf*2*pi;           % ��Ƶ��Ƶ��

HALF.fre_hc = HALF.f_rf * HALF.n_hc;  % г��ǻδʧгг��Ƶ��
HALF.wre_hc = HALF.fre_hc * 2*pi;     % г��ǻδʧгг���Ƶ��
HALF.qc = HALF.T0*HALF.I0/HALF.h/HALF.fillrate;         % bunch �����
HALF.wr_hc = HALF.w_rf*HALF.n_hc+HALF.fre_shift_hc*pi*2;% г��ǻʧгг���Ƶ��
HALF.angle_hc = HALF.wr_hc * HALF.Tb;     % г��ǻ���� norminal bucket center ��ת��

HALF.fre_hc2 = HALF.f_rf * HALF.n_hc2;  % г��ǻ2δʧгг��Ƶ��
HALF.wre_hc2 = HALF.fre_hc2 * 2*pi;     % г��ǻ2δʧгг���Ƶ��
HALF.wr_hc2 = HALF.w_rf*HALF.n_hc2+HALF.fre_shift_hc2*pi*2;% г��ǻ2ʧгг���Ƶ��
HALF.angle_hc2 = HALF.wr_hc2 * HALF.Tb;     % г��ǻ2���� norminal bucket center ��ת��

HALF.fre_mc = HALF.f_rf;                  % ��ǻδʧгг��Ƶ��
HALF.wre_mc = HALF.fre_mc * 2*pi;         % ��ǻδʧгг���Ƶ��
HALF.wr_mc = HALF.wre_mc + HALF.fre_shift_mc*pi*2;
HALF.angle_mc = HALF.wr_mc * HALF.Tb;     % ��ǻ���� norminal bucket center ��ת��

HALF.lambda_rads = HALF.T0 / HALF.tau_s;  % ����������   ��λ��Ȧ��
HALF.lambda_radz = HALF.T0 / HALF.tau_z;  % ���Ӽ�����   ��λ��Ȧ��

% Bunch form factor
% HALF.Fac_hc_init = Factor_calc(HALF.sigma_t0,HALF.fre_hc,HALF.Q_hc,HALF.wre_hc);
% HALF.Fac_hc2_init = Factor_calc(HALF.sigma_t0,HALF.fre_hc2,HALF.Q_hc2,HALF.wre_hc2);
% HALF.Fac_mc_init = Factor_calc(HALF.sigma_t0,HALF.fre_mc,HALF.Q_mc,HALF.wre_mc);
HALF.Fac_hc_init = F1;
HALF.Fac_hc2_init = F2;
HALF.Fac_mc_init = F0;


HALF.fais_nat = pi-asin(HALF.U0/HALF.V_mc);% natural synchrotron phase
HALF.det_angle_hc= atan(HALF.Q_hc*(HALF.wr_hc/HALF.wre_hc-HALF.wre_hc/HALF.wr_hc));
disp(['����г��ǻ���ز��ٶ�������������ʧг�Ƕȣ�',num2str(HALF.det_angle_hc/pi*180),' deg']);
HALF.fais_mc_whc = pi - asin((HALF.U0+2*HALF.I0*HALF.R_hc*cos(HALF.det_angle_hc)^2)/HALF.V_mc);
disp(['����г��ǻ���ز��ٶ�������������ͬ����λ��',num2str(HALF.fais_mc_whc),' rad  ',num2str(HALF.fais_mc_whc/pi*180),' deg']);
HALF.V_hc_load_0 = 2*HALF.Fac_hc_init*HALF.I0*HALF.R_hc*cos(HALF.det_angle_hc)*exp(1i*HALF.det_angle_hc);
disp(['����г��ǻ���ز��ٶ������������ĳ�ʼ���أ�',num2str(HALF.V_hc_load_0)]);

HALF.det_angle_hc2= atan(HALF.Q_hc2*(HALF.wr_hc2/HALF.wre_hc2-HALF.wre_hc2/HALF.wr_hc2));
disp(['����г��ǻ2���ز��ٶ�������������ʧг�Ƕȣ�',num2str(HALF.det_angle_hc2/pi*180),' deg']);
HALF.V_hc2_load_0 = 2*HALF.Fac_hc2_init*HALF.I0*HALF.R_hc2*cos(HALF.det_angle_hc2)*exp(1i*HALF.det_angle_hc2);
disp(['����г��ǻ2���ز��ٶ������������ĳ�ʼ���أ�',num2str(HALF.V_hc2_load_0)]);

HALF.det_angle_mc= atan(HALF.Q_mc*(HALF.wr_mc/HALF.wre_mc-HALF.wre_mc/HALF.wr_mc));
disp(['������ǻ���ز��ٶ�������������ʧг�Ƕȣ�',num2str(HALF.det_angle_mc/pi*180),' deg']);
HALF.V_mc_load_0 = 2*HALF.Fac_mc_init*HALF.I0*HALF.R_mc*cos(HALF.det_angle_mc)*exp(1i*HALF.det_angle_mc);
disp(['������ǻ���ز��ٶ������������ĳ�ʼ���أ�',num2str(HALF.V_mc_load_0)]);

HALF.Vmc_set = -HALF.V_mc*exp(-1i*(pi/2-HALF.fais_mc)); % MCǻѹ����λ��Ϊ����ֵ
HALF.Vg_mc_init = HALF.Vmc_set-HALF.V_mc_load_0;
HALF.Vg_mc_track_0 = HALF.Vg_mc_init; % initial generator voltage phasor
disp(['��ǻ�������ѹ��',num2str(HALF.Vg_mc_init)]);

HALF.Ig_mc_0 = HALF.Vg_mc_track_0/(2*HALF.R_mc*cos(HALF.det_angle_mc))*exp(-1i*HALF.det_angle_mc);
HALF.Vg_mc_0 = 2*pi*HALF.R_mc/HALF.Q_mc*HALF.Ig_mc_0;

HALF.drift_coef = HALF.alpha_c * HALF.T0 * (HALF.sigma_e0 / HALF.sigma_t0);
HALF.kick_coef  = 1 / (HALF.E0 * HALF.sigma_e0);

HALF.rfcoef1 = HALF.kick_coef * HALF.V_mc; % HALF.V_mc need be changed later

HALF.rfcoef2 = HALF.w_rf * HALF.sigma_t0;
% HALF.rfcoef2 = (HALF.w_rf-30e3*2*pi) * HALF.sigma_t0; % �ı���ǻ��г��Ƶ��

HALF.radampcoef = 2 * HALF.lambda_rads;
HALF.quanexcoef = 2 * sqrt(HALF.lambda_radz);

HALF.ploss = HALF.kick_coef * HALF.U0;

HALF.VbImagFactor_hc   = 1/sqrt(4*Q_hc^2-1);
HALF.rot_coef_hc       = sqrt(1-1/(2*Q_hc)^2);

HALF.VbImagFactor_hc2   = 1/sqrt(4*Q_hc2^2-1);
HALF.rot_coef_hc2       = sqrt(1-1/(2*Q_hc2)^2);

HALF.VbImagFactor_mc   = 1/sqrt(4*Q_mc^2-1);
HALF.rot_coef_mc       = sqrt(1-1/(2*Q_mc)^2);

% ���㸺�ؽ�
HALF.angle_V_mc = Vb_angle_calc(real(HALF.Vmc_set),imag(HALF.Vmc_set));
HALF.angle_V_ge = Vb_angle_calc(real(HALF.Vg_mc_init),imag(HALF.Vg_mc_init));

HALF.angle_load = HALF.fais_mc-(pi/2 - HALF.angle_V_ge-HALF.det_angle_mc);
disp(['��ǻ��ʼ���ؽ�����Ϊ��',num2str(HALF.angle_load/pi*180),'deg']);


% active HHC
HALF.Vhc_set = -HALF.V_hc*exp(-1i*(pi/2-HALF.fais_hc));
HALF.Vg_hc_init = HALF.Vhc_set-HALF.V_hc_load_0;
HALF.Vg_hc_track_0 = HALF.Vg_hc_init; % initial generator voltage phasor
disp(['г��ǻ�������ѹ��',num2str(HALF.Vg_hc_init)]);

HALF.Ig_hc_0 = HALF.Vg_hc_track_0/(2*HALF.R_hc*cos(HALF.det_angle_hc))*exp(-1i*HALF.det_angle_hc);
HALF.Vg_hc_0 = 2*pi*HALF.n_hc*HALF.R_hc/HALF.Q_hc*HALF.Ig_hc_0;

HALF.angle_V_hc = Vb_angle_calc(real(HALF.Vhc_set),imag(HALF.Vhc_set));
HALF.angle_V_ge_hc = Vb_angle_calc(real(HALF.Vg_hc_init),imag(HALF.Vg_hc_init));
% HALF.angle_V_ge_hc = angle(HALF.Vg_hc_init);

% HALF.angle_load_hc = pi/2+HALF.fais_hc-(HALF.angle_V_ge_hc-HALF.det_angle_hc);
HALF.angle_load_hc = HALF.fais_hc-(pi/2 - HALF.angle_V_ge_hc-HALF.det_angle_hc);
disp(['г��ǻ��ʼ���ؽ�����Ϊ��',num2str(HALF.angle_load_hc/pi*180),'deg']);

% active HHC2
HALF.Vhc2_set = -HALF.V_hc2*exp(-1i*(pi/2-HALF.fais_hc2));
HALF.Vg_hc2_init = HALF.Vhc2_set-HALF.V_hc2_load_0;
HALF.Vg_hc2_track_0 = HALF.Vg_hc2_init; % initial generator voltage phasor
disp(['г��ǻ2�������ѹ��',num2str(HALF.Vg_hc2_init)]);

HALF.Ig_hc2_0 = HALF.Vg_hc2_track_0/(2*HALF.R_hc2*cos(HALF.det_angle_hc2))*exp(-1i*HALF.det_angle_hc2);
HALF.Vg_hc2_0 = 2*pi*HALF.n_hc2*HALF.R_hc2/HALF.Q_hc2*HALF.Ig_hc2_0;

HALF.angle_V_hc2 = Vb_angle_calc(real(HALF.Vhc2_set),imag(HALF.Vhc2_set));
HALF.angle_V_ge_hc2 = Vb_angle_calc(real(HALF.Vg_hc2_init),imag(HALF.Vg_hc2_init));
% HALF.angle_V_ge_hc2 = angle(HALF.Vg_hc_init);

%HALF.angle_load_hc2 = pi/2+HALF.fais_hc-(HALF.angle_V_ge-HALF.det_angle_hc);
HALF.angle_load_hc2 = HALF.fais_hc2-(pi/2 - HALF.angle_V_ge_hc2-HALF.det_angle_hc2);
disp(['г��ǻ2��ʼ���ؽ�����Ϊ��',num2str(HALF.angle_load_hc2/pi*180),'deg']);

end