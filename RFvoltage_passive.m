function [V_load_kick] = RFvoltage_passive(HALF,exp_ang_coef,TbAng_coef,pattern,V_load_0_real,V_load_0_imag,Q,i,id)
    exp_angle  = exp(exp_ang_coef * Q);
    exp_angle_sum= gather(sum(exp_angle));       % 耗时 0.007s sum()函数较慢    

    if id == 1
        V_load_cpu = double(exp_angle_sum.*HALF.Vb_mc(pattern==1));   
        [V_load_real,V_load_imag,V_load_0_real,V_load_0_imag] = VoltageLoadCalc(real(TbAng_coef),...
        imag(TbAng_coef),V_load_0_real,V_load_0_imag,real(V_load_cpu),imag(V_load_cpu),pattern);
    else 
        if id == 2
            V_load_cpu = double(exp_angle_sum.*HALF.Vb_hc(pattern==1));    
            [V_load_real,V_load_imag,V_load_0_real,V_load_0_imag] = VoltageLoadCalc(real(TbAng_coef),...
            imag(TbAng_coef),V_load_0_real,V_load_0_imag,real(V_load_cpu),imag(V_load_cpu),pattern);
        else
            V_load_cpu = double(exp_angle_sum.*HALF.Vb_hc2(pattern==1));    
            [V_load_real,V_load_imag,V_load_0_real,V_load_0_imag] = VoltageLoadCalc(real(TbAng_coef),...
            imag(TbAng_coef),V_load_0_real,V_load_0_imag,real(V_load_cpu),imag(V_load_cpu),pattern);
        end
    end 
    
    V_load_cpu = (V_load_real + 1i * V_load_imag)*HALF.kick_coef;   % 约化V_load;
    V_load     = gpuArray(single(V_load_cpu));    
    % intrabunch kick    - V_load_kick    real part
    V_load_kick = V_load./exp_angle;

    if mod(i,1000)==0
        figure(id)
        subplot(2,1,1)
        plot(abs(V_load_real + 1i * V_load_imag)/1e3);
        title('Harmonic Cavity: Beam Loading');ylabel('Amplitude [kV]');
        subplot(2,1,2)
        plot(angle(V_load_real + 1i * V_load_imag)-pi/2);ylabel('Phase [deg]');xlabel('Bucket ID');
        subplot(3,1,3)
        plot(V_load_real/1e3);ylabel('real(V) [kV]');xlabel('Bucket ID');
    end

end

