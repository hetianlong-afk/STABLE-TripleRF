function [Response_Vb]=Response_matrix(h,omega_r,R_hc,Q_hc,T0)
Response_Vb = zeros(h);
V0 = omega_r*R_hc/Q_hc;
V1 = V0/(1-exp((1i-1/(2*Q_hc))*omega_r*T0));
Vb0 = zeros(1,h);
Vb0(1) = V1*exp((1i-1/(2*Q_hc))*omega_r*T0);
for i = 2 : h
    Vb0(i) = V1*exp((1i-1/(2*Q_hc))*omega_r*(i-1)*T0/h);
end
for i = 1 : h
    Response_Vb(i,i:end) = Vb0(1:end-i+1);
     if i>1
        Response_Vb(i,1:i-1)= Vb0(end-i+2:end);
     end
end

%for i = 1 : h
%    for ii = 1 : h
%        if i < ii
%            Response_Vb(i,ii) = V1*exp((1i-1/2/Q_hc)*omega_r*(ii-i)*T0/h);
%        else
%            Response_Vb(i,ii) = V1*exp((1i-1/2/Q_hc)*omega_r*(ii-i+h)*T0/h);
%        end
%    end
%end


%Vb = zeros(1,h);


%for i = 1 : h
%    for ii = 1 : h
%        if pattern(ii) == 1
%            Vb(i) = Vb(i)+Response_Vb(ii,i)*T0*I0/Bun_num*F(ii);
%        end
%    end
%end

%Bun_index=1:1:h;
%plot(Bun_index,abs(Vb)*1e-3);ylabel('Induced voltage [kV]');xlabel('Bunch index');

end
