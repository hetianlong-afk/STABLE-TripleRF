function [V_kick] = RFvoltage_ideal(HALF,exp_ang_coef,pattern,Q,id)
    if id == 1
        V = HALF.V_mc;
        fais = HALF.fais_mc;
    else 
        if id == 2
        V = HALF.V_hc;
        fais = HALF.fais_hc;
        else
        V = HALF.V_hc2;
        fais = HALF.fais_hc2;
        end
    end

    exp_angle  = exp(exp_ang_coef * Q);
    Vc = -V*exp(-1i*(pi/2-fais))*ones(1,HALF.h);
    V_kick = gpuArray(single(Vc(pattern==1)*HALF.kick_coef))./exp_angle;
end

