%%  统计作密度分布图  Dq = 0.4;
Dq = 0.2;
Q_min = min(Q);
tau_min = gather(Q_min)*HALF.sigma_t0;
Q_new = round((Q - Q_min) *(1/Dq));
binnum=max(max(Q_new))+1; binnum=gather(binnum);
bin_num_q=sum(BinNumCalZ(binnum,Q_new));
bin_num_q=reshape(bin_num_q,binnum,Bun_num);

% figure(4);
% bin_i = [1:5:100,320:5:420,640];
% bin_i = [1:4:25,130:4:154,260:4:284];
% bin_i = 1;
bin_i = [1:HALF.h];
% bin_i = [1:160:HALF.h];
% bin_i = [1:10:100];
colorset =[1 0 0;0 1 0;0 0 1;0 0 0.5;1 0.5 0.5; 0.5 1 0.5;0.5 0.5 1;1 0 1;0 1 1;1 1 0];
tau = zeros(binnum,length(bin_i));
for i =1:length(bin_i)
    bin_range = (1:binnum)*Dq*HALF.sigma_t0+tau_min(bin_i(i));
    % plot((bin_range'-Dq*HALF.sigma_t0)*1e12,bin_num_q(:,bin_i(i))/1e4,'-','Linewidth',2);hold on; % rs
    tau(:,i) = (bin_range'-Dq*HALF.sigma_t0);
end
% ylabel('norm. density [a.u.]');xlabel('\tau [ps]'); 
% xlim([-250,250]);
% set(gca,'FontName','Times New Roman','FontSize',14);

delta_fais = record_Q_mean(end,:)*HALF.sigma_t0*HALF.w_rf;

% Fac1 = zeros(1,HALF.h);
% Fac2 = zeros(1,HALF.h);
% for i =1:length(bin_i)
%     Fac1(i) = gather(sum(bin_num_q(:,bin_i(i))/1e4.*exp(-1i*tau(:,i)*HALF.wr_hc)));
%     Fac2(i) = gather(sum(bin_num_q(:,bin_i(i))/1e4.*exp(-1i*tau(:,i)*HALF.wr_hc2)));
% end

F1 = 0.826435778739579 + 0.003158864487639i;
F2 = 0.567832151901265 + 0.003414962774284i;

t = 1 : HALF.h;
F1_real = 1e-6*sin(t/HALF.h*2*pi)+real(F1);
F2_real = 1e-6*sin(t/HALF.h*2*pi)+real(F2);
% F1_real = real(F1)*ones(1,HALF.h);
% F2_real = real(F2)*ones(1,HALF.h);
F1_imag = 1i*1e-6*sin(t/HALF.h*2*pi)+1i*imag(F1);
F2_imag = 1i*1e-6*sin(t/HALF.h*2*pi)+1i*imag(F2);
% F1_imag = 1i*imag(F1)*ones(1,HALF.h);
% F2_imag = 1i*imag(F2)*ones(1,HALF.h);
Fac1 = F1_real+F1_imag;
Fac2 = F2_real+F2_imag;

% 
figure(5)
subplot(2,1,1)
plot(real(Fac1));hold on;
subplot(2,1,2)
plot(imag(Fac1));hold on;
figure(6)
subplot(2,1,1)
plot(real(Fac2));hold on;
subplot(2,1,2)
plot(imag(Fac2));hold on;


Fac1_new = repmat(Fac1.',1,HALF.h);
Fac2_new = repmat(Fac2.',1,HALF.h);

Response_hc = Response_matrix(HALF.h,HALF.wr_hc,HALF.R_hc,HALF.Q_hc,HALF.T0);  % SSRM
Response_hc2 = Response_matrix(HALF.h,HALF.wr_hc2,HALF.R_hc2,HALF.Q_hc2,HALF.T0);

% lambda_hc = eig(Response_hc*HALF.qc(1)*1e4);
% lambda_hc2 = eig(Response_hc2*HALF.qc(1)*1e4);
% figure(10)
% subplot(2,1,1)
% plot(abs(lambda_hc));hold on;
% subplot(2,1,2)
% plot(abs(angle(lambda_hc)));hold on;
% figure(11)
% subplot(2,1,1)
% plot(abs(lambda_hc2));hold on;
% subplot(2,1,2)
% plot(abs(angle(lambda_hc2)));hold on;
% V_load_hc = sum(Response_hc.*exp(-1i*HALF.n_hc*(delta_fais'-delta_fais)).*Fac1_new*HALF.qc(1)*1e4);  % 稳态负载
% V_load_hc2 = sum(Response_hc2.*exp(-1i*HALF.n_hc2*(delta_fais'-delta_fais)).*Fac2_new*HALF.qc(1)*1e4);  % 稳态负载
V_load_hc = sum(Response_hc.*Fac1_new*HALF.qc(1)*1e4);  % 稳态负载
V_load_hc2 = sum(Response_hc2.*Fac2_new*HALF.qc(1)*1e4);  % 稳态负载




figure(1)
subplot(3,1,1)
plot(abs(V_load_hc));hold on;
subplot(3,1,2)
plot(angle(V_load_hc));hold on;
subplot(3,1,3)
plot(real(V_load_hc));hold on;
figure(2)
subplot(3,1,1)
plot(abs(V_load_hc2));hold on;
subplot(3,1,2)
plot(angle(V_load_hc2));hold on;
subplot(3,1,3)
plot(real(V_load_hc2));hold on;

% for i = 1 : binnum
%     lambda0(i) = 1/sqrt(2*pi)/HALF.sigma_t0*exp(-((i-round(binnum/2))*1e-12)^2/2/(HALF.sigma_t0^2))/1e12;
%     t(i) = (i-round(binnum/2))*1e-12;
% end
% den_sum = sum(lambda0);
% mean_tau = sum(t.*lambda0)/den_sum;
% sigma_tau = sqrt(sum(lambda0.*(t-mean_tau).^2)/den_sum);
% 
% Fac0_1 = sum(lambda0.*exp(-1i*t*HALF.wr_hc));
% Fac0_2 = sum(lambda0.*exp(-1i*t*HALF.wr_hc2));
% 
% figure(3)
% plot(t*1e12,lambda0);hold on;

