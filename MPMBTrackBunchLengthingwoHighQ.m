% 多粒子多束团纵向追踪
% 研究谐波腔致束团拉伸
% GPU 加速
% 作者： 何天龙
% 时间： 20200706
% time:  20200910 (modified: w/o high Q approximation)
% time : 20210213 (modified: include the BL of main cavity)
% time : 20211121 (modified: include HOMs)
% time : 20220212 (modified: include realistic PI feedback for MC)
% time : 20240115 (modified: include realistic PI feedback for HHC)
clc;clear;
%% beam parameters
% HALF 参数
cspeed = 299792458; 
sigma_t0 = 10e-12;      % s initial rms bunch length  （用于计算上归一化）
sigma_e0 = 7.44e-4;     %  rms energy spread
alpha_c  = 9.4e-5;      %  momentum compaction factor
tau_s    = 14e-3;       %  radiation damping time
tau_z    = 14e-3;       %  radiation damping time
I0     = 350e-3;        %  beam current
E0     = 2.2e9;         %  beam energy
U0     = 400e3;         %  energy loss per turn
V_mc   = 1.2e6;         %  main cavity voltage
h      = 800;           % harmonic number
n_hc1  = 3;             % harmonic order of HHC1
n_hc2  = 5;             % harmonic order of HHC2
C      = 479.86;        % Circumference of the ring
% Q_hc   = 1e5;           % HHC loaded quality factor
% R_hc   = Q_hc*40;       % HHC loaded shunt impedance   
Q_hc2  = 1e5;           % HHC loaded quality factor
R_hc2  = Q_hc2*2;      % HHC loaded shunt impedance   
Q_mc   = 1e5;           % shunt impedance of main cavity
R_mc   = Q_mc*45;       % quality factor of main cavity
%% RF cavity parameter settings for different longitudinal bunch equilibrium distribution
% obtained by analytical formula

fais_mc = 2.74028221659595;    %%1.2MV flat potential condition
fais_hc1 = -0.140514110229679;
fais_hc2 = 3.05692778011039;
k1 = 0.464856515347598;
k2 = 0.0923859093693406;

F0 = 0.979626801531162 + 0.001275780712274i; %%1.2MV flat potential condition
F1 = 0.826435778739579 + 0.003158864487639i;
F2 = 0.567832151901265 + 0.003414962774284i;


V_hc = k1 * V_mc;   %ideal voltage amplitude of lower order HC
V_hc2 = k2 * V_mc;  %ideal voltage amplitude of higher order HC
%% RF cavity detuning setting (Minimum generator power)

T0 = C/cspeed;
f_rf = h/T0;

psi_hc = angle(F1)-fais_hc1+pi/2;                            %detuning angle of lower order HC
R_hc = k1*V_mc/(abs(-2*I0*abs(F1)*cos(psi_hc)*exp(-1i*psi_hc)));
Q_hc = R_hc/15;                                              
fre_shift = -n_hc1*f_rf*tan(psi_hc)/2/Q_hc;                  

psi_hc2 = atan(abs(F2)*2*I0*R_hc2*cos(fais_hc2)/k2/V_mc);    %detuning angle of higer order HC
fre_shift2 = tan(psi_hc2)*h*cspeed*n_hc2/C/2/Q_hc2;          

psi_mc = atan(abs(F0)*2*I0*R_mc*cos(fais_mc)/V_mc);          %detuning angle of MC
fre_shift_mc = tan(psi_mc)*h*cspeed/C/2/Q_mc;
%% Working mode setting

% fill pattern
pattern = zeros(1,h);
pattern(1:1:h)=1;    %Uniform filling
fillrate = length(find(pattern==1))/h;

HALF = machine(C,I0,U0,E0,tau_s,tau_z,sigma_t0,sigma_e0,alpha_c,h,...
    V_mc,V_hc,V_hc2,fais_mc,fais_hc1,fais_hc2,n_hc1,n_hc2,R_hc,Q_hc...
    ,R_hc2,Q_hc2,fillrate,fre_shift,fre_shift2,Q_mc,R_mc,fre_shift_mc,F1,F2,F0);
HALF.ShortRange_on = 0; % 0 - neglecting short range effect, 1 - considering.

HALF.MC_mode  = 1; % 1 - ideal RF cavity, 2 - passive RF cavity, 3 - active RF cavity.
HALF.HC_mode  = 3; % 1 - ideal RF cavity, 2 - passive RF cavity, 3 - active RF cavity.
HALF.HC2_mode  = 3; % 1 - ideal RF cavity, 2 - passive RF cavity, 3 - active RF cavity.

PI_Set;   % add PI 
%% bunch generation
Par_num = 1e4; Bun_num = length(find(pattern==1));

%% charge pattern

% HALF
charge = ones(1,h).*pattern; 

q1 = TruncatedGaussian(1, [-3,3], [Par_num,1]);
p1 = TruncatedGaussian(1, [-3,3], [Par_num,1]);
q  = repmat(q1,1,Bun_num);
p  = repmat(p1,1,Bun_num);%
%CPU to GPU
Q=gpuArray(single(q)); P=gpuArray(single(p));     % single type

index_add = 1:Bun_num;
index_add = gpuArray(single(index_add-1));

Dq = 0.05;
%% wake data (intrabunch motion)
tau_q = (0:Dq:200)'*sigma_t0;
Wake_inter_hc = -HALF.wr_hc *  R_hc /Q_hc*exp(-tau_q*HALF.wr_hc/2/Q_hc) .*(...
    cos(tau_q*HALF.wr_hc*HALF.rot_coef_hc)-HALF.VbImagFactor_hc*sin(tau_q*HALF.wr_hc*HALF.rot_coef_hc));
Wake_inter_hc(1) = Wake_inter_hc(1)/2;

Wake_inter_mc = -HALF.wr_mc *  R_mc /Q_mc*exp(-tau_q*HALF.wr_mc/2/Q_mc) .*(...
    cos(tau_q*HALF.wr_mc*HALF.rot_coef_mc)-HALF.VbImagFactor_mc*sin(tau_q*HALF.wr_mc*HALF.rot_coef_mc));
Wake_inter_mc(1) = Wake_inter_mc(1)/2;

Wake_inter = Wake_inter_hc + Wake_inter_mc;

bin_tau = Dq*sigma_t0;

Wake_inter    = gpuArray(single(Wake_inter));
%% start tracking Track_num = 1e3
% charge per macro-particle   : HALF.qc
% 不等电荷量填充时，每个束团的宏粒子电荷量不等，注意区别

HALF.qc   = charge.*pattern * HALF.qc / Par_num;              %由单个元素变为一行矩阵
% induced voltage per macro-particle  : HALF.V_b
HALF.Vb_hc  = HALF.qc * HALF.wr_hc * HALF.R_hc / HALF.Q_hc; 
HALF.Vb_hc2 = HALF.qc * HALF.wr_hc2 * HALF.R_hc2 / HALF.Q_hc2; 
HALF.Vb_mc  = HALF.qc * HALF.wr_mc * HALF.R_mc / HALF.Q_mc; 

% HALF.V_b  = HALF.qc * HALF.w_r * HALF.R_hc / HALF.Q_hc *(1+1i*HALF.VbImagFactor); 
% initial loaded voltage

V_hc_load_0_real = real(HALF.V_hc_load_0);
V_hc_load_0_imag = imag(HALF.V_hc_load_0);
V_hc_load_0      = V_hc_load_0_real+1i*V_hc_load_0_imag;
% V_hc_load_0      =0;

V_hc2_load_0_real = real(HALF.V_hc2_load_0);
V_hc2_load_0_imag = imag(HALF.V_hc2_load_0);
V_hc2_load_0      = V_hc2_load_0_real+1i*V_hc2_load_0_imag;

V_mc_load_0_real = real(HALF.V_mc_load_0);
V_mc_load_0_imag = imag(HALF.V_mc_load_0);
V_mc_load_0      = V_mc_load_0_real+1i*V_mc_load_0_imag;
% V_mc_load_0=0; 


rot_decay_coef_hc = 1i * HALF.rot_coef_hc - 1 / (2 * HALF.Q_hc);  % 旋转 衰减项
TbAng_coef_hc     = exp(rot_decay_coef_hc * HALF.angle_hc);       % 注意符号正负
exp_ang_coef_hc   = -rot_decay_coef_hc * HALF.wr_hc * sigma_t0;
exp_ang_coef_hc   = gpuArray(single(exp_ang_coef_hc));

rot_decay_coef_hc2 = 1i * HALF.rot_coef_hc2 - 1 / (2 * HALF.Q_hc2);  % 旋转 衰减项
TbAng_coef_hc2     = exp(rot_decay_coef_hc2 * HALF.angle_hc2);       % 注意符号正负
exp_ang_coef_hc2   = -rot_decay_coef_hc2 * HALF.wr_hc2 * sigma_t0;
exp_ang_coef_hc2   = gpuArray(single(exp_ang_coef_hc2));

rot_decay_coef_mc = 1i * HALF.rot_coef_mc - 1 / (2 * HALF.Q_mc);  % 旋转 衰减项
TbAng_coef_mc     = exp(rot_decay_coef_mc * HALF.angle_mc);       % 注意符号正负
exp_ang_coef_mc   = -rot_decay_coef_mc * HALF.wr_mc * sigma_t0;
exp_ang_coef_mc   = gpuArray(single(exp_ang_coef_mc));

% Track number.
Track_num  = 30e4;

wake_kick_coef = HALF.qc * HALF.kick_coef;

% 发射机电压矢量替代之前的Vrf矢量
Vg_mc = abs(HALF.Vg_mc_init);
[Vg_angle]=round(Vb_angle_calc(real(HALF.Vg_mc_init),imag(HALF.Vg_mc_init))*1e12)/1e12;
HALF.Vg_mc_track = HALF.Vg_mc_init;
HALF.rfcoef1_track     = HALF.rfcoef1 / HALF.V_mc * Vg_mc;
fai_s_track            = pi/2-Vg_angle;                    % 发射机电压矢量的同步相位

% record parameters 
Recor_step = 100;
HALF.Recor_step=Recor_step;
Recor_num  = Track_num / Recor_step;
% Vg_mc_track_record = zeros(1,Track_num*h/10);
Vb_mc_track_record = zeros(1,Track_num);
Vb_hc_track_record = zeros(1,Track_num);
Vb_hc2_track_record = zeros(1,Track_num);
record_Q_mean = zeros(Recor_num,Bun_num);record_Q_std = zeros(Recor_num,Bun_num);
record_P_mean = zeros(Recor_num,Bun_num);record_P_std = zeros(Recor_num,Bun_num);
record_V_load = zeros(Recor_num,Bun_num);record_V_load2 = zeros(Recor_num,Bun_num);
%%
gd = gpuDevice(); tic;
j=1;
record_th = 0;

Ig_mc_track_tnum = 0;
Ig_hc_track_tnum = 0;
Ig_hc2_track_tnum = 0;
for i =1:Track_num
%%
    % drift
    Q = Q + P * HALF.drift_coef;   
    Q_min = min(Q);    
    Q_new = round((Q - Q_min) *(1/Dq));        
%% radiation damping and quantum excitation term   
    rad_quan_kick = -HALF.radampcoef * P + HALF.quanexcoef *...
        gpuArray.randn(Par_num,Bun_num,'single');
    P = P + rad_quan_kick - HALF.ploss;
%% Main cavity id = 1

    if HALF.MC_mode == 1
        [V_mc_kick] = RFvoltage_ideal(HALF,exp_ang_coef_mc,pattern,Q,1);
    else
        if HALF.MC_mode == 2
            [V_mc_kick] = RFvoltage_passive(HALF,exp_ang_coef_mc,TbAng_coef_mc,pattern,V_mc_load_0_real,V_mc_load_0_imag,Q,i,1);
        else
            [Vc_mc,V_load_mc,V_mc_kick,Vb_mc_track_record,PIMC,HALF] = RFvoltage_active(HALF,PIMC,exp_ang_coef_mc,TbAng_coef_mc,Vb_mc_track_record,pattern,Q,i,1);
        end
    end

    P = P - real(V_mc_kick) + imag(V_mc_kick) * HALF.VbImagFactor_mc;
% % _________________________________________________________________________   
%% Harmonic cavity  id = 2   
    if HALF.HC_mode == 1
        [V_hc_kick] = RFvoltage_ideal(HALF,exp_ang_coef_hc,pattern,Q,2);
    else
        if HALF.HC_mode == 2
            [V_hc_kick] = RFvoltage_passive(HALF,exp_ang_coef_hc,TbAng_coef_hc,pattern,V_hc_load_0_real,V_hc_load_0_imag,Q,i,2);
        else
            [Vc_hc,V_load_hc,V_hc_kick,Vb_hc_track_record,PIHC,HALF] = RFvoltage_active(HALF,PIHC,exp_ang_coef_hc,TbAng_coef_hc,Vb_hc_track_record,pattern,Q,i,2);
        end
    end
    
    P = P - real(V_hc_kick) + imag(V_hc_kick) * HALF.VbImagFactor_hc;
% _________________________________________________________________________  
%% Harmonic cavity2  id = 3  
    if HALF.HC2_mode == 1
        [V_hc2_kick] = RFvoltage_ideal(HALF,exp_ang_coef_hc2,pattern,Q,3);
    else
        if HALF.HC2_mode == 2
            [V_hc2_kick] = RFvoltage_passive(HALF,exp_ang_coef_hc2,TbAng_coef_hc2,pattern,V_hc2_load_0_real,V_hc2_load_0_imag,Q,i,3);
        else
            [Vc_hc2,V_load_hc2,V_hc2_kick,Vb_hc2_track_record,PIHC2,HALF] = RFvoltage_active(HALF,PIHC2,exp_ang_coef_hc2,TbAng_coef_hc2,Vb_hc2_track_record,pattern,Q,i,3);
        end
    end

    P = P - real(V_hc2_kick) + imag(V_hc2_kick) * HALF.VbImagFactor_hc2;
% _________________________________________________________________________ 
%%  short range wakefield
    % count bins
    if HALF.ShortRange_on ==1                 % modified in 2022/11/14
        binnum=max(max(Q_new))+1; binnum=gather(binnum);
        bin_num_q=sum(BinNumCalZ(binnum,Q_new));   % double type
        bin_num_q=reshape(bin_num_q,binnum,Bun_num);    
        kick_conv = conv2(bin_num_q,Wake_inter(1:binnum));    
        kick_conv = kick_conv(1:binnum,:) .* wake_kick_coef(pattern==1)*min(i/5000,1);
    % wake kick    + wake_kick
        Q_new = Q_new + (1 + index_add * binnum);  % modified in 2020/09/22
        wake_kick = kick_conv(Q_new);
    else
        wake_kick = 0;
    end
    % radiation damping and quantum excitation term  + wake_kick
  
    P = P + wake_kick;
%%    
    if mod(i,2000)==0
        Centroid_std=std(record_Q_mean(record_th,:))*HALF.sigma_t0*1e12;
        disp(['tracking turn = ',num2str(i),'; Centroid_std = ',num2str(Centroid_std),'ps']);
        toc;
    end
%    output data
    if mod(i,Recor_step)==0
        record_th = record_th +1;
        record_Q_mean(record_th,:)=gather(mean(Q));
        record_Q_std(record_th,:)=gather(std(Q));
        record_P_mean(record_th,:)=gather(mean(P));
        record_P_std(record_th,:)=gather(std(P));
    end
end
wait(gd);
toc;
%%
filename=[num2str(Par_num/1e4),'_',num2str(R_hc2/Q_hc2),'_',num2str(R_hc/Q_hc),'_','PI.m',num2str(PIHC.m),'_',num2str(HALF.tau_s),'_',num2str(Track_num),'_',num2str(Recor_step),'_',num2str(I0*1e3),'.mat'];
save(filename,'record_Q_mean','record_Q_std','record_P_mean','record_P_std','Q','Track_num','HALF','PIMC','PIHC','PIHC2','Bun_num');
%%
figure(1);
Recor_step=HALF.Recor_step;

Nturns = (1:Track_num/Recor_step)*Recor_step;
for i=1:50:HALF.h
    subplot(2,2,1)
    plot(Nturns,record_Q_mean(:,i)*HALF.sigma_t0*1e12); hold on;
    subplot(2,2,2)
    plot(Nturns,record_Q_std(:,i)*HALF.sigma_t0*1e12); hold on;
    subplot(2,2,3)

    plot(Nturns,record_P_mean(:,i)*HALF.sigma_e0); hold on;
    subplot(2,2,4)
    plot(Nturns,record_P_std(:,i)*HALF.sigma_e0); hold on; 
%     pause(0.2);
end
% Track_num=50000;
subplot(2,2,1);ylabel('<\tau>  [ps]');xlabel('turns');xlim([1,Track_num]);grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
subplot(2,2,2);ylabel('\sigma_{\tau}  [ps]');xlabel('turns');xlim([1,Track_num]);grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
subplot(2,2,3);ylabel('<\delta> ');xlabel('turns');xlim([1,Track_num]);grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
subplot(2,2,4);ylabel('\sigma_{\delta} ');xlabel('turns');xlim([1,Track_num]);grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
%%
% 统计沿着束团 长度分布，中心分布  1:100:2000
figure(2);
t_plot = 1:HALF.h;
Q_mean_plot = zeros(1,HALF.h); 
Q_std_plot = zeros(1,HALF.h);
Q_mean_plot(pattern==1) = record_Q_mean(end,:);
Q_mean_plot(pattern==0) = NaN;
Q_std_plot(pattern==1) = record_Q_std(end,:);
Q_std_plot(pattern==0) = NaN;
subplot(1,2,1);plot(Q_mean_plot*HALF.sigma_t0*1e12,'.','MarkerSize',7);hold on;
ylabel('<\tau>  [ps]');xlabel('bunch number');
set(gca,'FontName','Times New Roman','FontSize',10.8);
subplot(1,2,2);plot(Q_std_plot*HALF.sigma_t0*1e12,'.','MarkerSize',7);hold on;
ylabel('\sigma_{\tau}  [ps]');xlabel('bunch number');
set(gca,'FontName','Times New Roman','FontSize',10.8);
%%  统计作密度分布图  Dq = 0.4;
Dq = 0.2;
Q_min = min(Q);
tau_min = gather(Q_min)*HALF.sigma_t0;
Q_new = round((Q - Q_min) *(1/Dq));
binnum=max(max(Q_new))+1; binnum=gather(binnum);
bin_num_q=sum(BinNumCalZ(binnum,Q_new));
bin_num_q=reshape(bin_num_q,binnum,Bun_num);

delta_tau = 1.6667*1e-12;
figure(4);
bin_i = [400];
% colorset =[1 0 0;0 1 0;0 0 1;0 0 0.5;1 0.5 0.5; 0.5 1 0.5;0.5 0.5 1;1 0 1;0 1 1;1 1 0];
for i =1:length(bin_i)
    bin_range = (1:binnum)*Dq*HALF.sigma_t0+tau_min(bin_i(i));
    plot((bin_range'-Dq*HALF.sigma_t0)*1e12,bin_num_q(:,bin_i(i))/Par_num,'.','Markersize',10);hold on; % rs
end

ylabel('norm. density [a.u.]');xlabel('\tau [ps]'); 
set(gca,'FontName','Times New Roman','FontSize',12);
%% 画出Vg电压

figure(656)
% subplot(2,1,1)
plot(abs(Vb_hc_track_record)/1e3/mean(abs(Vb_hc_track_record(1:50000))/1e3));title('Harmonic Cavity');hold on;
% subplot(2,1,2)
plot(angle(Vb_hc_track_record)/pi*180/mean(angle(Vb_hc_track_record(1:50000))/pi*180));hold on;
% plot(angle(Vb_hc_track_record)/pi*180);
plot(Nturns,record_P_mean(:,1)/7+1);hold on;
legend('Amplitude','Phase','<\delta>');
xlabel('Turns');ylabel('Norm.Amp. [a.u.]');xlim([0,10e4]);
%%
%% 主腔发射机功率  PI.Ig_track 
% 发射机电流
figure(667)
subplot(2,1,1);
plot(abs(PIMC.Ig_track));hold on
plot(abs(PIHC.Ig_track));ylabel('Ig amplitude [A]');
legend('MC','HC');
subplot(2,1,2);
plot(angle(PIMC.Ig_track)/pi*180);hold on;
plot(angle(PIHC.Ig_track)/pi*180);ylabel('Ig phase [deg]');