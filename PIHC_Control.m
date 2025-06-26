function [Vc_hc,Vg_hc_track,Vg_hc_track_0,V_load,V_hc_load_0,PIHC]=PIHC_Control(PIHC,Vrf_ideal,Vg_hc_track_0,...
    V_hc_load_0,V_load_cpu,TbAng_coef_hc,pattern)
% realistic PI feedback control
% PI.m 每m个buckets的腔压平均值  For HALF, PI.m=5120, PI.d=800
% PI.d 延迟d个buckets输出反馈量
% Ig=Ig0+DIg
% i is turn number
% HALF.Vhc_set <- Vrf_ideal
% HALF.Vg_hc_track_0
% Vc_mc : total voltage phasor = generator + beam driven
% V_load : beam loading voltage phasor
n          =length(pattern);
V_load     =zeros(1,n);
Vg_hc_track=zeros(1,n);
Vc_hc      =zeros(1,n);
j=1;
for ii =1:n
    PIHC.piIndex = PIHC.piIndex+1; % 第 piIndex 个buckets,依据其值施加反馈
%% PI feedback   PI.d > 0 is an integer 
    % if PIHC.piIndex>PIHC.d && mod(PIHC.piIndex-PIHC.d,PIHC.m)==0  % N*m+d 时施加反馈
    %    PIHC.Vg0 = 2*pi*PIHC.n_hc*PIHC.RoverQ*PIHC.Ig_FBid(1);        % impluse phasor single pass
	%    PIHC.Ig_FBid(1)=[];
    % end
% 注释此段表示不改变Vg0
%% RF phase modulation
% to test
%     if PIHC.piIndex > 5e4*n
%         t_offset = RFModulationTimeOffest(PIHC.RFModulation_mf,...
%             PIHC.RFModulation_wm,PIHC.RFModulation_wf,PIHC.RFModulation_Tb,PIHC.piIndex);       
% %         PIHC.Vg0 = 2*pi*PIHC.RoverQ*(PIHC.Ig0*exp(1i*t_offset*PIHC.RFModulation_wf)-PIHC.DIg);
%         PIHC.Vg0 = 2*pi*PIHC.RoverQ*(PIHC.Ig*exp(-1i*t_offset*PIHC.RFModulation_wf));
%     end
%%    
    Vg_hc_track_0  = Vg_hc_track_0 * TbAng_coef_hc + PIHC.Vg0;
    Vg_hc_track(ii)= Vg_hc_track_0;
    if pattern(ii)==1
        V_load(ii) = V_hc_load_0 * TbAng_coef_hc;
        V_hc_load_0 = V_load(ii) + V_load_cpu(j);
        j=j+1;
    else
        V_hc_load_0 = V_hc_load_0 * TbAng_coef_hc;
        V_load(ii) = V_hc_load_0;
    end
    Vc_hc(ii) = Vg_hc_track_0 + V_load(ii);    %Total voltage phasor
%% PI feedback   
    piTrackIndex = mod(PIHC.piIndex,PIHC.m);
    if piTrackIndex==0 %&& PIHC.piIndex > 1e0*n  % N*m 时计算腔压平均矢量,并输送给PI计算反馈输出量
        piTrackIndex = PIHC.m;
        PIHC.Vrf_hc_track(piTrackIndex) = Vc_hc(ii);%
        Vrf_hc_mean = mean(PIHC.Vrf_hc_track);
        DI = (Vrf_ideal-Vrf_hc_mean)/PIHC.RL; % RL 负载阻抗   此处前后关系
        PIHC.Integral=PIHC.Integral+PIHC.KI*DI;   % DI error signal
        PIHC.DIg = PIHC.KP*DI+PIHC.Integral;
        PIHC.Ig = PIHC.Ig0 + PIHC.DIg;            % 发射机电流，注意此处为'+'
        PIHC.Ig_track = [PIHC.Ig_track,PIHC.Ig];
		PIHC.Ig_FBid  = [PIHC.Ig_FBid,PIHC.Ig];
        PIHC.Vg0 = 2*pi*PIHC.n_hc*PIHC.RoverQ*PIHC.Ig_FBid(1);  % impluse phasor single pass
	    PIHC.Ig_FBid(1)=[];
    end   
    if piTrackIndex==0
        piTrackIndex = PIHC.m;
    end    
    PIHC.Vrf_hc_track(piTrackIndex) = Vc_hc(ii);%
end
end