function [Vc_mc,Vg_mc_track,Vg_mc_track_0,V_load,V_mc_load_0,PI]=PI_Control_copy64(PI,Vrf_ideal,Vg_mc_track_0,...
    V_mc_load_0,V_load_cpu,TbAng_coef_mc,pattern)
% realistic PI feedback control : CASE 2
% PI.m ÿm��buckets��ǻѹƽ��ֵ  For HALF, PI.m=5120, PI.d=800
% PI.d �ӳ�d��buckets��������� FIR �������64�ݣ����θ�PI 
% Ig=Ig0+DIg
% i is turn number
% HALF.Vrf_ideal
% HALF.Vg_mc_track_0
% Vc_mc : total voltage phasor = generator + beam driven
% V_load : beam loading voltage phasor
n          =length(pattern);
V_load     =zeros(1,n);
Vg_mc_track=zeros(1,n);
Vc_mc      =zeros(1,n);
j=1;
for ii =1:n
    PI.piIndex = PI.piIndex+1; % �� piIndex ��buckets,������ֵʩ�ӷ���
    
    if PI.piIndex < 8e4*n   % ���� or �ջ�
        PIFBON=1;
    else
        PIFBON=0;
    end
%% PI feedback   PI.d > 0 is an integer 
    if PI.piIndex>PI.m && mod(PI.piIndex-PI.d,80)==0 && PIFBON==1  % N*m+d ʱʩ�ӷ���
       PI.Integral=PI.Integral+PI.KI*PI.DI;   % DI error signal
       PI.DIg = PI.KP*PI.DI+PI.Integral;
       PI.Ig = PI.Ig0 - PI.DIg;
       PI.Vg0 = 2*pi*PI.RoverQ*PI.Ig;         % impluse phasor single pass
%        disp(['piIndex = ',num2str(PI.piIndex),' ',num2str(PI.Vg0)]);%test
    end
% ע�ʹ˶α�ʾ���ı�Vg0
%% RF phase modulation
% to test
    if PI.piIndex > 5e4*n
        t_offset = RFModulationTimeOffest(PI.RFModulation_mf,...
            PI.RFModulation_wm,PI.RFModulation_wf,PI.RFModulation_Tb,PI.piIndex);       
        PI.Vg0 = 2*pi*PI.RoverQ*(PI.Ig*exp(-1i*t_offset*PI.RFModulation_wf));
    end
%%    
    Vg_mc_track_0  = Vg_mc_track_0 * TbAng_coef_mc + PI.Vg0;
    Vg_mc_track(ii)= Vg_mc_track_0;
    if pattern(ii)==1
        V_load(ii) = V_mc_load_0 * TbAng_coef_mc;
        V_mc_load_0 = V_load(ii) + V_load_cpu(j);
        j=j+1;
    else
        V_mc_load_0 = V_mc_load_0 * TbAng_coef_mc;
        V_load(ii) = V_mc_load_0;
    end
    Vc_mc(ii) = Vg_mc_track_0 + V_load(ii);     %Total voltage phasor
%% PI feedback   
    piTrackIndex = mod(PI.piIndex,PI.m);
    if piTrackIndex==0 && PI.piIndex > 1e4*n  % N*m ʱ����ǻѹƽ��ʸ��,�����͸�PI���㷴�������
        piTrackIndex = PI.m;
        PI.Vrf_mc_track(piTrackIndex) = Vc_mc(ii);%
        Vrf_mc_mean = mean(PI.Vrf_mc_track);
        PI.DI = (Vrf_mc_mean-Vrf_ideal)/PI.RL; % RL �����迹
%         disp(['PI.DIg = ',num2str(PI.DIg)]); %test
    end   
    if piTrackIndex==0
        piTrackIndex = PI.m;
    end    
    PI.Vrf_mc_track(piTrackIndex) = Vc_mc(ii);%    
end
end