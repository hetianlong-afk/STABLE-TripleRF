function [Vc,V_load,V_kick,Vb_track_record,PI,HALF] = RFvoltage_active(HALF,PI,exp_ang_coef,TbAng_coef,Vb_track_record,pattern,Q,i,id)
    exp_angle  = exp(exp_ang_coef * Q);
    exp_angle_sum= gather(sum(exp_angle));       % 耗时 0.007s sum()函数较慢   

    if id == 1
        V_load_cpu = double(exp_angle_sum.*HALF.Vb_mc(pattern==1));  
        [Vc,Vg_mc_track,HALF.Vg_mc_track_0,V_load,HALF.V_mc_load_0,PI]=PIMC_Control(PI,HALF.Vmc_set,HALF.Vg_mc_track_0,...
        HALF.V_mc_load_0,V_load_cpu,TbAng_coef,pattern); % 
    else 
        if id == 2
            V_load_cpu = double(exp_angle_sum.*HALF.Vb_hc(pattern==1));   
            [Vc,Vg_hc_track,HALF.Vg_hc_track_0,V_load,HALF.V_hc_load_0,PI]=PIHC_Control(PI,HALF.Vhc_set,HALF.Vg_hc_track_0,...
            HALF.V_hc_load_0,V_load_cpu,TbAng_coef,pattern); % 
        else
            V_load_cpu = double(exp_angle_sum.*HALF.Vb_hc2(pattern==1));   
            [Vc,Vg_hc2_track,HALF.Vg_hc2_track_0,V_load,HALF.V_hc2_load_0,PI]=PIHC2_Control(PI,HALF.Vhc2_set,HALF.Vg_hc2_track_0,...
            HALF.V_hc2_load_0,V_load_cpu,TbAng_coef,pattern); % 
        end
    end 

    Vb_track_record(i)=mean(V_load);
	Ig_hc_track_num = length(PI.Ig_track);
    if Ig_hc_track_num>5000
        PI.Ig_track = [];
    end
%  谐波腔腔压矢量图示
    if mod(i,1000)==0
        figure(id)
        subplot(3,1,1)
        plot(abs(V_load)/1e3);
        title('Harmonic Cavity: Beam Loading');ylabel('Amplitude [kV]');
        subplot(3,1,2)
        plot(angle(V_load)-pi/2);ylabel('Phase [deg]');xlabel('Bucket ID');
        subplot(3,1,3)
        plot(real(V_load)/1e3);ylabel('real(V) [kV]');xlabel('Bucket ID');
        figure(id*11)
        subplot(2,1,1)
        plot(abs(Vc)/1e3);
        title('Harmonic Cavity: Total Voltage');ylabel('Amplitude [kV]');
        subplot(2,1,2)
        plot(angle(Vc)-pi/2);ylabel('Phase [deg]');xlabel('Bucket ID');
		figure(id*111);
        subplot(2,1,1);plot(abs(PI.Ig_track));title('HC Generator current');ylabel('Amplitude');
        subplot(2,1,2);plot(angle(PI.Ig_track));ylabel('Phase'); 
    end 

	V_kick = gpuArray(single(Vc(pattern==1)*HALF.kick_coef))./exp_angle;
end

