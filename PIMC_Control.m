
function [Vc_mc,Vg_mc_track,Vg_mc_track_0,V_load,V_mc_load_0,PIMC]=PIMC_Control(PIMC,Vrf_ideal,Vg_mc_track_0,...
    V_mc_load_0,V_load_cpu,TbAng_coef_mc,pattern)
% realistic PI feedback control
% PI.m ÿm��buckets��ǻѹƽ��ֵ  For HALF, PI.m=5120, PI.d=800
% PI.d �ӳ�d��buckets���������
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
    PIMC.piIndex = PIMC.piIndex+1; % �� piIndex ��buckets,������ֵʩ�ӷ���
%% PI feedback   PI.d > 0 is an integer 
    % if PIMC.piIndex>PIMC.d && mod(PIMC.piIndex-PIMC.d,PIMC.m)==0  % N*m+d ʱʩ�ӷ���
    %    PIMC.Vg0 = 2*pi*PIMC.RoverQ*PIMC.Ig_FBid(1);        % impluse phasor single pass
	%    PIMC.Ig_FBid(1)=[];
    % end
% ע�ʹ˶α�ʾ���ı�Vg0
%% RF phase modulation
% to test
%     if PIMC.piIndex > 5e4*n
%         t_offset = RFModulationTimeOffest(PIMC.RFModulation_mf,...
%             PIMC.RFModulation_wm,PIMC.RFModulation_wf,PIMC.RFModulation_Tb,PIMC.piIndex);       
% %         PIMC.Vg0 = 2*pi*PIMC.RoverQ*(PIMC.Ig0*exp(1i*t_offset*PIMC.RFModulation_wf)-PIMC.DIg);
%         PIMC.Vg0 = 2*pi*PIMC.RoverQ*(PIMC.Ig*exp(-1i*t_offset*PIMC.RFModulation_wf));
%     end
%%    
    Vg_mc_track_0  = Vg_mc_track_0 * TbAng_coef_mc + PIMC.Vg0;
    Vg_mc_track(ii)= Vg_mc_track_0;
    if pattern(ii)==1
        V_load(ii) = V_mc_load_0 * TbAng_coef_mc;
        V_mc_load_0 = V_load(ii) + V_load_cpu(j);
        j=j+1;
    else
        V_mc_load_0 = V_mc_load_0 * TbAng_coef_mc;
        V_load(ii) = V_mc_load_0;
    end
    Vc_mc(ii) = Vg_mc_track_0 + V_load(ii);    %Total voltage phasor
%% PI feedback   
    piTrackIndex = mod(PIMC.piIndex,PIMC.m);
    if piTrackIndex==0 %&& PIMC.piIndex > 1e0*n  % N*m ʱ����ǻѹƽ��ʸ��,�����͸�PI���㷴�������
        piTrackIndex = PIMC.m;
        PIMC.Vrf_mc_track(piTrackIndex) = Vc_mc(ii);%
        Vrf_mc_mean = mean(PIMC.Vrf_mc_track);
        DI = (Vrf_ideal-Vrf_mc_mean)/PIMC.RL; % RL �����迹   �˴�ǰ���ϵ
        PIMC.Integral=PIMC.Integral+PIMC.KI*DI;   % DI error signal
        PIMC.DIg = PIMC.KP*DI+PIMC.Integral;
        PIMC.Ig = PIMC.Ig0 + PIMC.DIg;            % �����������ע��˴�Ϊ'+'
        PIMC.Ig_track = [PIMC.Ig_track,PIMC.Ig];
		PIMC.Ig_FBid  = [PIMC.Ig_FBid,PIMC.Ig];
        PIMC.Vg0 = 2*pi*PIMC.RoverQ*PIMC.Ig_FBid(1); % impluse phasor single pass
	    PIMC.Ig_FBid(1)=[];
    end   
    if piTrackIndex==0
        piTrackIndex = PIMC.m;
    end    
    PIMC.Vrf_mc_track(piTrackIndex) = Vc_mc(ii);%
end
end