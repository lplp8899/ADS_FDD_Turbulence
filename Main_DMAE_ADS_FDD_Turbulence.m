%% ----------------------------------------------------------------------
% This is the code for the following journal paper:
%
%      P. Lu, E. van Kampen, C. C. de Visser, Q. P. Chu
%      Air Data Sensor Fault Detection and Diagnosis in the Presence of Atmospheric Turbulence: 
%      Theory and Experimental Validation with Real Flight Data.
%      IEEE transactions on Control Systems Technology, 2020.
%
%   This program shows the results of Section IV of the paper.
%
%   Author: Peng Lu, the University of Hong Kong
%
%   email:  lupeng@hku.hk
%
%   released in July 2020
%
%   Please cite the above paper if you find this code useful
%----------------------------------------------------------------------

%%
tic;

clear all
close all
% clc

% addpath('/Data');
addpath(fullfile(pwd, 'Data'));

% choose whether you want to have the Q adaptation or not
disp('0. do NOT use Q adaption     1. use Q adaption')
Q_flag = input('Please choose: ');

% choose different turbulence scenarios, the bigger number, the more
% intensive the turbulence is
disp('0. no tur    1.tur 1    2. tur 2    3. tur 3   4.tur 4')
case_flag=input('please choose one tur:  ');

% simultaneous or multiple faults
disp('0.multiple fault     1.simultanous faults  ')
simul_flag=input('please choose one situation:  ');


%% ----------------------- five tubulence scenarios ----------------------
% case 0 is no turbulence, which is not included in the paper (Table 1)
if case_flag == 0
    load new_sim_Ma_08_H_5000_tur0.mat
elseif case_flag == 1
    load new_sim_Ma_08_H_5000_tur1.mat
elseif case_flag == 2
    load new_sim_Ma_08_H_5000_tur2.mat 
elseif case_flag == 3
    load new_sim_Ma_08_H_5000_tur3.mat 
elseif case_flag == 4
    load new_sim_Ma_08_H_5000_tur4.mat 
end


%% sampling time
delta_t = 0.01;

% r2d = 180/pi;


%% get true states from the data
Vtas   =  Out_nonoise(:,1)';
alpha  =  Out_nonoise(:,2)';
beta   =  Out_nonoise(:,3)';

pb     =  Out_nonoise(:,4)';
qb     =  Out_nonoise(:,5)';
rb     =  Out_nonoise(:,6)';

phi    =  Out_nonoise(:,9)';
theta  =  Out_nonoise(:,8)';
psi    =  Out_nonoise(:,7)';

xe     =  Out_nonoise(:,10)';
ye     =  Out_nonoise(:,11)';
ze     =  Out_nonoise(:,12)'; % opposite sign of he
u_b    =  Out_nonoise(:,13)';
v_b    =  Out_nonoise(:,14)';
w_b    =  Out_nonoise(:,15)';

u_n    =  Out_nonoise(:,16)';
v_n    =  Out_nonoise(:,17)';
w_n    =  Out_nonoise(:,18)'; % the opposite sign of hdot


% IMU measurements with noise
A_x    =  squeeze(Ax_m)';
A_y    =  squeeze(Ay_m)';
A_z    =  squeeze(Az_m)';

% turbulence
ug_b     = Tur_b(:,1)';
vg_b     = Tur_b(:,2)';
wg_b     = Tur_b(:,3)';


% lower the dimension
Out_noise=squeeze(Out_noise)';
% get measurenments from the data
p_m         =  Out_noise(:,4)';
q_m         =  Out_noise(:,5)';
r_m         =  Out_noise(:,6)';

Ax_m        =  squeeze(Ax_m)'; % due to the simulation
Ay_m        =  squeeze(Ay_m)';
Az_m        =  squeeze(Az_m)';

Vt_m        =  Out_noise(:,1)';
alpha_m     =  Out_noise(:,2)';
beta_m      =  Out_noise(:,3)';

phi_m       =  Out_noise(:,9)';
theta_m     =  Out_noise(:,8)';
psi_m       =  Out_noise(:,7)';
ze_m        =  Out_noise(:,12)'; % opposite sign of he_m
% ze          =  Out_noise(:,17)';
xe_m        =  Out_noise(:,10)';
ye_m        =  Out_noise(:,11)';

%
u_n_m       =  Out_noise(:,16)';
v_n_m       =  Out_noise(:,17)';
w_n_m       =  Out_noise(:,18)'; % Vz= -hedot, Vz and wn are in the same direction


%% initialisation for the filter
x_ob_0 = zeros(12,1);
% x_ob_0 = [xe_m(1); ye_m(1); ze_m(1); u_b(1); v_b(1); w_b(1); phi_m(1); theta_m(1); psi_m(1); zeros(6,1)];
x_ob_0 = [Vt_m(1); alpha_m(1); beta_m(1); phi_m(1); theta_m(1); psi_m(1); ];

P_x0 =  1e2*eye(6);

       
% IMU measurements as the input
u  =  [Ax_m; Ay_m; Az_m; p_m; q_m; r_m];
[dim_u,gen] = size(u);


% true state vector from simulation
x_real = [ Vtas; alpha; beta; phi; theta; psi;];

% real measurement vector
z_real = [ Vt_m;  alpha_m;  beta_m;  phi_m; theta_m;  psi_m;  ];


%% ------------------------ inject sensors faults start ------------------------
z_real_p=z_real;
% 
% % order of the faults
n_fault1=1;
n_fault2=2;
n_fault3=3;
%% ----------------------------   multiple faults, but not simultaneous   --------------------------
if simul_flag==0
fy=0.0001;
% epsu_s=0.3;
faults_sensor=zeros(3,gen);
faults_sensor(n_fault1,1000:2000)=faults_sensor(n_fault1,1000:2000)+ 4*5/4*1;
faults_sensor(n_fault2,3000:4000)=faults_sensor(n_fault2,3000:4000)+ 2/57.3*1;
faults_sensor(n_fault3,5000:6000)=faults_sensor(n_fault3,5000:6000)+ 2/57.3*1*1;

z_real(1,1:gen)=z_real(1,1:gen)+faults_sensor(1,:);
z_real(2,1:gen)=z_real(2,1:gen)+faults_sensor(2,:);
z_real(3,1:gen)=z_real(3,1:gen)+faults_sensor(3,:);

%-------------------------   Ideal PDF for plotting -----------------------------
PI=zeros(2,gen);
PI(1,1:999)=1;PI(1,2001:2999)=1;PI(1,4001:4999)=1;PI(1,6001:8000)=1;
PI(2,1000:2000)=1;PI(2,3000:4000)=1;PI(2,5000:6000)=1;
%-------------------------   Ideal PDF for plotting -----------------------------
end


%% ----------------------------------   simultaneous faults   ----------------------------
if simul_flag == 1
fy=0.0001;
faults_sensor=zeros(3,gen);
%
%------------ in the presence of normal turbulence ------------
faults_sensor(n_fault1,1000:2000)=faults_sensor(n_fault1,1000:2000)+ 2*3/2;
% faults_sensor(n_fault2,1000:2000)=faults_sensor(n_fault2,1000:2000)+ 0.5/57.3;
faults_sensor(n_fault2,1000:2000)=faults_sensor(n_fault2,1000:2000)+ 1/57.3*2;
faults_sensor(n_fault3,1000:2000)=faults_sensor(n_fault3,1000:2000)+ 1/57.3*2;
%
faults_sensor(n_fault1,3000:4000)=faults_sensor(n_fault1,3000:4000)- 2*3/2;
faults_sensor(n_fault2,3000:4000)=faults_sensor(n_fault2,3000:4000)- 1/57.3*2;
faults_sensor(n_fault3,3000:4000)=faults_sensor(n_fault3,3000:4000)- 1/57.3*2;
%------------ in the presence of big turbulence ------------

z_real(1,1:gen)=z_real(1,1:gen)+faults_sensor(1,:);
z_real(2,1:gen)=z_real(2,1:gen)+faults_sensor(2,:);
z_real(3,1:gen)=z_real(3,1:gen)+faults_sensor(3,:);

%-------------------------   Ideal PDF for plotting -----------------------------
PI=zeros(2,gen);
PI(1,1:999)=1;PI(1,2001:2999)=1;PI(1,4001:8000)=1;
PI(2,1000:2000)=1;PI(2,3000:4000)=1;%PI(4,5000:6000)=1;
%-------------------------   Ideal PDF for plotting -----------------------------
end
%------------------------ inject sensors faults end   ------------------------



% acceleration of gravity
% g = y(:,28)';
g = 9.81*ones(1,gen);

nameda_factor = 1;

%-------------------------------------------------------------------------

%% ------------------------------  DMAE initialization parameters  ------------------------------
[dim_sys,~] = size(x_ob_0);
[dim_out,~] = size(z_real);
%
x_ob = zeros(dim_sys,2*dim_sys+1);
z_ob = zeros(dim_out,2*dim_sys+1);
x_ob_mean = zeros(dim_sys,gen);
z_ob_mean = zeros(dim_out,gen);
inno = zeros(dim_out,gen);
norm_inno = zeros(1,gen);
error_x = zeros(6,gen);
norm_error_x = zeros(dim_sys,gen);
x_ob_filter = zeros(dim_sys,gen);
z_ob_filter = zeros(dim_out,gen);
P_ob_filter = zeros(dim_sys,dim_sys,gen);
residual = zeros(dim_out,gen);
norm_residual = zeros(1,gen);

Qk = diag([1e-4,1e-4,1e-4,3e-8,3e-8,3e-8]);
Rk = diag([1e-2,3e-6,3e-6,3e-8,3e-8,3e-8]); 

% Qk = diag([(10^(-1))^2,(10^(-1))^2,(10^(-1))^2,3e-8,3e-8,3e-8]);

% initial for af
mean_x0=x_ob_0;
mean_x0_ff=mean_x0;
P_x0_ff=1e0*eye(6);
mean_x0_af=[mean_x0; 0; 0; 0;];
P_x0_af=1e0*eye(9);
%
[dim_sys_af,~]=size(mean_x0_af);

inno_ff=zeros(dim_out,gen);
residual_ff = zeros(dim_out,gen);
error_x_ff=zeros(6,gen);
x_ob_filter_ff=zeros(dim_sys,gen);
z_ob_filter_ff =zeros(dim_out,gen);
Pzk_ff=zeros(dim_out,dim_out,gen);
P_ob_filter_ff=zeros(dim_sys,dim_sys,gen);

inno_af=zeros(dim_out,gen);
residual_af = zeros(dim_out,gen);
error_x_af=zeros(6,gen);
x_ob_filter_af=zeros(dim_sys_af,gen);
z_ob_filter_af =zeros(dim_out,gen);
Pzk_af=zeros(dim_out,dim_out,gen);
P_ob_filter_af=zeros(dim_sys_af,dim_sys_af,gen);

%
x_ob_filter_pdf = zeros(dim_sys,gen);
x_ob_filter_pdf2 = zeros(dim_sys,gen);

% --------------------- est ---------------------
Pk_ff = zeros(dim_sys,dim_sys,gen);
Pk_af = zeros(dim_sys_af,dim_sys_af,gen);
S_inno_ff=zeros(dim_out,dim_out,gen);
C_inno_ff=zeros(dim_out,dim_out,gen);
S_inno_af=zeros(dim_out,dim_out,gen);
C_inno_af=zeros(dim_out,dim_out,gen);

Qd_est_ff = zeros(dim_sys,dim_sys,gen); 
Qd_est_af = zeros(dim_sys_af,dim_sys_af,gen); 
Qfk = zeros(3,3);

Kk_ff = zeros(dim_sys,dim_out,gen);
Kk_af = zeros(dim_sys_af,dim_out,gen);
% --------------------- est ---------------------

pd_ff=zeros(1,gen);
pd_af=zeros(1,gen);

p_t=zeros(1,gen);
p_ff=zeros(1,gen);
p_af=zeros(1,gen);

p_k=zeros(2,gen);

liho_ff=zeros(1,gen);
liho_af=zeros(1,gen);

% pdf min & max
I_min=zeros(1,gen);
I_max=zeros(1,gen);

% the initial PDF
p0_0=0.9;
p1_0=0.1;

K_flag=[];
%
Kf=[];
Knf=[];

f_ob_sq = zeros(3,3,gen);

% parameters for UKF
ks=0;% k = 3 - n if n < 3
as=0.8;% as<1
% as=0.5;% as<1
bs=2;
namedas=as.^2*(dim_sys+ks)-dim_sys;
w_m0=namedas/(dim_sys+namedas);% weight coeffiency for mean value
w_c0=namedas/(dim_sys+namedas)+(1-as.^2+bs);% weight coeffiency for mean value
w_mi=1/(2*(dim_sys+namedas));
w_ci=w_mi;
% column vector
w_m=[w_m0; w_mi*ones(2*dim_sys,1);];
w_c=[w_c0; w_ci*ones(2*dim_sys,1);];
mid_m=(eye(2*dim_sys+1)-w_m*ones(1,2*dim_sys+1));
W_mc=mid_m*diag(w_c')*mid_m';

w_m0_af=namedas/(dim_sys_af+namedas);% weight coeffiency for mean value
w_c0_af=namedas/(dim_sys_af+namedas)+(1-as.^2+bs);% weight coeffiency for mean value
w_mi_af=1/(2*(dim_sys_af+namedas));
w_ci_af=w_mi_af;
w_m_af=[w_m0_af; w_mi_af*ones(2*dim_sys_af,1);];


%% iteration of the DMAE approach
for k=1:gen/(gen/6000)
    
    % run the fault free filter
    [inno_ff(:,k),residual_ff(:,k),error_x_ff(:,k),Kk_ff(:,:,k),Pk_ff(:,:,k),Pzk_ff(:,:,k),x_ob_filter_ff(:,k),P_ob_filter_ff(:,:,k),z_ob_filter_ff(:,k)]...
        = UKF_fault_free(mean_x0_ff,P_x0_ff,x_real,z_real,Qk,Rk,u,g,delta_t,namedas,w_c0,w_ci,w_m,dim_sys,k,Qd_est_ff(:,:,k));

    % run the augmented fault filter
    [inno_af(:,k),residual_af(:,k),error_x_af(:,k),Kk_af(:,:,k),Pk_af(:,:,k),Pzk_af(:,:,k),x_ob_filter_af(:,k),P_ob_filter_af(:,:,k),z_ob_filter_af(:,k)]...
        = UKF_augmented_fault(mean_x0_af,P_x0_af,x_real,z_real,Qk,Rk,u,g,delta_t,namedas,w_c0_af,w_ci_af,w_m_af,dim_sys_af,k,Qd_est_af(:,:,k));


    %% ------------------------------ PDF update -----------------------
    % likelihood function
    liho_ff(:,k)=(inno_ff(:,k)'/(2*Pzk_ff(:,:,k))*inno_ff(:,k));
    liho_af(:,k)=(inno_af(:,k)'/(2*Pzk_af(:,:,k))*inno_af(:,k));

    pd_ff(:,k)=exp(-liho_ff(:,k))/(det(2*pi*Pzk_ff(:,:,k))^(1/2));
    pd_af(:,k)=exp(-liho_af(:,k))/(det(2*pi*Pzk_af(:,:,k))^(1/2));

    % iteration
    pd0=pd_ff(:,k)*p0_0;
    p1=pd_af(:,k)*p1_0;

    % sum of the PDF
    p_t(:,k)=pd0+p1;
    p_ff(:,k)=pd0/p_t(:,k);
    p_af(:,k)=p1/p_t(:,k);
  
    %% --------------  initial values for next iteration     ---------------
    mean_x0_ff=x_ob_filter_ff(:,k);
    P_x0_ff=P_ob_filter_ff(:,:,k);

    mean_x0_af=x_ob_filter_af(:,k);
    P_x0_af=P_ob_filter_af(:,:,k);
    
  
    % mark the index for max & min pdf
    p_k(:,k)=[p_ff(:,k);p_af(:,k)];
    [min_p,I_min(k)]=min(p_k(:,k));
    [max_p,I_max(k)]=max(p_k(:,k));
    
    flag_f=0;  

    %% ------------------------- start of the selective reinitialization -------------------
    % ------------  no faults:   reinitialise aug filter using non-aug filter  ------------
    % if p_ff(:,k) > p_af(:,k)
    if I_max(k)==1 %&& sum(kf)<2

        Knf=[Knf;k];
        mean_x0_af(1:dim_sys,:)=mean_x0_ff;  P_x0_af(1:dim_sys,1:dim_sys)=P_x0_ff;
        mean_x0_af(dim_sys+1:dim_sys_af)=1e-3*ones(3,1);
        P_x0_af(1:dim_sys,dim_sys+1:dim_sys_af)=zeros(dim_sys,3);
        P_x0_af(dim_sys+1:dim_sys_af,1:dim_sys)=zeros(3,dim_sys);       
        P_x0_af(dim_sys+1:dim_sys_af,dim_sys+1:dim_sys_af)=1e2*eye(3);
    end

    % --------------  faults:   reinitialise non-aug filter using aug filter    ------------
    if k>1 && I_max(k)==2 %&& flag_f~=1%&& (I_max(k-1)==2) && (I_max(k-5)==2) 
        Kf=[Kf;k];
        mean_x0_ff=mean_x0_af(1:dim_sys,:);
        P_x0_ff=P_x0_af(1:dim_sys,1:dim_sys);
    end
    % ------------------------- end of the selective reinitialization -------------------
   
    
    %% min, prevent lock out
    p_ff(:,k)=max(p_ff(:,k),0.001);
    p_af(:,k)=max(p_af(:,k),0.001);

    % max, prevent lock out
    p_ff(:,k)=min(p_ff(:,k),0.999);
    p_af(:,k)=min(p_af(:,k),0.999);
    
    % the next 
    p0_0=p_ff(:,k);
    p1_0=p_af(:,k);
    
    
    %----------------------- fusion of estimates ------------------------
    x_ob_filter_pdf(:,k)=p_ff(:,k)*x_ob_filter_ff(:,k)+p_af(:,k)*x_ob_filter_af(1:dim_sys,k);
   
    % choose the model with the highest PDF, do not fuse all the models
    [p_max,I]=max([p_ff(:,k),p_af(:,k)]);
    x_com=[x_ob_filter_ff(:,k),x_ob_filter_af(1:dim_sys,k)];
    x_ob_filter_pdf2(:,k)=x_com(:,I);


    %% -------------------------------  Q adaptation starts -------------------------------
    %--------------------- innovation covariance est -------------------
    S_inno_af(:,:,k) = inno_af(:,k)*inno_af(:,k)';
 
    N_win = 10;% width of moving window
    N_win = 40;% width of moving window
    
    % sum of the innovation square
    if N_win == 1
        C_inno_af(:,:,k) = S_inno_af(:,:,k);
    end
    if N_win > 1
        if k == 1
            C_inno_af(:,:,k) = S_inno_af(:,:,k);
        elseif k >1 && k < N_win
            C_inno_af(:,:,k) = sum(S_inno_af(:,:,1:k),3)/(k-1);
        elseif k > N_win
            C_inno_af(:,:,k) = sum(S_inno_af(:,:,k-N_win+1:k),3)/(N_win-1);
        end 
    end
       
    % -------------------------- Q estimation approximation -------------
    % use Fault Model to est
    G = G_Ja(x_ob_filter_af(:,k));

%     Qd_0 = C_inno_af(:,:,k) -  G*Qk*G'*delta_t - Rk;
    Fk = [eye(3); zeros(3,3)];
    Qd_0 = C_inno_af(:,:,k) -  G*Qk*G'*delta_t - Fk*Qfk*Fk'*delta_t - Rk;

    % guarantee positive for the covariance
    Qd_est = diag([ max(Qd_0(1,1),0), max(Qd_0(2,2),0), max(Qd_0(3,3),0),...
                    max(Qd_0(4,4),0), max(Qd_0(5,5),0), max(Qd_0(6,6),0)]);
    % Qd for fault free filter and augmente fault filter
    Qd_est_ff(:,:,k+1) = diag([ Qd_est(1,1), Qd_est(2,2)*1, Qd_est(3,3)*1, 0, 0, 0]);
    Qd_est_af(:,:,k+1) = diag([ Qd_est(1,1), Qd_est(2,2)*1, Qd_est(3,3)*1, 0, 0, 0, zeros(1,3) ]);

    if (Q_flag == 0) % if you do not use the Q adaptation, you can use fixed Q
        % Qd is non-zero, the paper used this one
        % note that this is ad hoc.
        Qd_est_ff(:,:,k+1) = diag([ 0, 0, 0, 0, 0, 0]);
        Qd_est_af(:,:,k+1) = diag([ 0, 0, 0, 0, 0, 0, 1e-4, (1e-2/57.3)^2, (1e-2/57.3)^2]);
        % Qd is zero
%         Qd_est_ff(:,:,k+1) = diag([ 0, 0, 0, 0, 0, 0]);
%         Qd_est_af(:,:,k+1) = diag([ 0, 0, 0, 0, 0, 0, 0, 0, 0]);

    end
    %-------------------------------  Q adaptation ends -------------------------------
    

%% covariance
f_ob = x_ob_filter_af(7:9,k);
f_ob_sq(:,:,k) = f_ob*f_ob'; % covariance

end


%% --------------  pdf filtering stage, not very necessary  ------------------%
p_ff_filter = p_ff;
p_af_filter = p_af;
for i = 1:k
    if p_ff_filter(i)>0.5
        p_ff_filter(i) = 0.999;
    end
    if p_ff_filter(i)<0.5
        p_ff_filter(i) = 0.001;
    end
    if p_af_filter(i)>0.5
        p_af_filter(i) = 0.999;
    end
    if p_af_filter(i)<0.5
        p_af_filter(i) = 0.001;
    end
end
% --------------  pdf filtering stage  ------------------%

%% --------------- plots

Time=delta_t*(1:k);


% turbulence speeds
figure;hold on;
subplot(311);plot(Time,Tur_b(1:k,1),'b','linewidth',1.2);ylabel('u_w (m/s)','fontsize',14); grid; set(gca,'fontsize',14);
% title('Turbulence')
subplot(312);plot(Time,Tur_b(1:k,2),'b','linewidth',1.2);ylabel('v_w (m/s)','fontsize',14); grid; set(gca,'fontsize',14);
subplot(313);plot(Time,Tur_b(1:k,3),'b','linewidth',1.2);ylabel('w_w (m/s)','fontsize',14); grid; set(gca,'fontsize',14);
h1=axes('position',[0.50 0.0001 0.0001 0.0001],'fontsize',14); title('time (s)','fontsize',14)

% Qd estimation
figure; axis off;
subplot(311);hold on; plot(Time(2:k),squeeze(Qd_est_af(1,1,2:k)),'b','linewidth',1.2); ylabel('Q_{d,11}^0','fontsize',14);grid;set(gca,'fontsize',14);
% title('Qd estimation, af')
subplot(312);hold on; plot(Time(2:k),squeeze(Qd_est_af(2,2,2:k)),'b','linewidth',1.2);  ylabel('Q_{d,22}^0','fontsize',14);grid;set(gca,'fontsize',14);
subplot(313);hold on; plot(Time(2:k),squeeze(Qd_est_af(3,3,2:k)),'b','linewidth',1.2); ylabel('Q_{d,33}^0','fontsize',14);grid;set(gca,'fontsize',14);
h1=axes('position',[0.50 0.0001 0.0001 0.0001],'fontsize',14); title('time (s)','fontsize',14)
% set(h2,'Box','off')


% model probabilities
figure;hold on;axis off;
subplot(211);hold on; plot(Time,PI(1,1:k),'r','linewidth',2); plot(Time,p_ff(1:k),'b--','linewidth',2); ylabel('p_{nf}','fontsize',14);set(gca,'fontsize',14);
% box off; set(gca,'position',[0.1 0.76 0.85 0.15],'fontsize',14);
subplot(212);hold on; plot(Time,PI(2,1:k),'r','linewidth',2); plot(Time,p_af(1:k),'b--','linewidth',2); ylabel('p_{af}','fontsize',14);set(gca,'fontsize',14);
% h1=axes('position',[0.09 0.05 0.88 0.85],'fontsize',14);
% axis off;
% title(h1,(situ(num)));
h=legend('True','DMAE'); %set(h,'color','none','edgecolor','white');
xlabel('time (s)');

% filtered model probabilities, almost the same with unfiltered. 
figure;hold on;axis off;
subplot(211);hold on; plot(Time,PI(1,1:k),'r','linewidth',2); plot(Time,p_ff_filter(1:k),'b--','linewidth',2); ylabel('p_{nf}','fontsize',14);set(gca,'fontsize',14);
% box off; set(gca,'position',[0.1 0.76 0.85 0.15],'fontsize',14);
subplot(212);hold on; plot(Time,PI(2,1:k),'r','linewidth',2); plot(Time,p_af_filter(1:k),'b--','linewidth',2); ylabel('p_{af}','fontsize',14);set(gca,'fontsize',14);
% h1=axes('position',[0.09 0.05 0.88 0.85],'fontsize',14);
% axis off;
% title(h1,(situ(num)));
h=legend('True','DMAE'); %set(h,'color','none','edgecolor','white');
xlabel('time (s)');


% estimates of the faults
figure;
subplot(311); hold on; plot(Time,faults_sensor(1,1:k),'r','linewidth',2); plot(Time,x_ob_filter_af(7,1:k),'b--','linewidth',2);  ylabel('f_{V} (m/s)','fontsize',14);grid; set(gca,'fontsize',14);% set(gca,'xlim',[0 delta_t*k],'ylim',[-2 12],'fontsize',14);
% title('Fault estimation')
subplot(312); hold on; plot(Time,faults_sensor(2,1:k),'r','linewidth',2); plot(Time,x_ob_filter_af(8,1:k),'b--','linewidth',2); ylabel('f_{\alpha} (rad)','fontsize',14);grid; set(gca,'fontsize',14);  %set(gca,'xlim',[0 delta_t*k],'ylim',[-0.05 0.05],'fontsize',14);
subplot(313); hold on; plot(Time,faults_sensor(3,1:k),'r','linewidth',2); plot(Time,x_ob_filter_af(9,1:k),'b--','linewidth',2); ylabel('f_{\beta} (rad)','fontsize',14);grid; set(gca,'fontsize',14); % set(gca,'xlim',[0 delta_t*k],'ylim',[-0.05 0.05],'fontsize',14);
h2=legend('True','DMAE'); set(h2,'color','white','edgecolor','black');
% h1=axes('position',[0.10 0.05 0.84 0.86],'fontsize',14);axis off;title(h1,'Fault estimation')
h1=axes('position',[0.50 0.0001 0.0001 0.0001],'fontsize',14); title('time (s)','fontsize',14)
% set(h2,'Box','off')

% probabiliti-weighted fault estimation
figure;
subplot(311); hold on; plot(Time,faults_sensor(1,1:k),'r','linewidth',2); plot(Time,p_af(1:k).*x_ob_filter_af(7,1:k),'b--','linewidth',2);  ylabel('f_{V} (m/s)','fontsize',14);grid; set(gca,'fontsize',14);% set(gca,'xlim',[0 delta_t*k],'ylim',[-0.2 1.2],'fontsize',14);
% title('Probability-weighted Fault estimation')
subplot(312); hold on; plot(Time,faults_sensor(2,1:k),'r','linewidth',2); plot(Time,p_af(1:k).*x_ob_filter_af(8,1:k),'b--','linewidth',2); ylabel('f_{\alpha} (rad)','fontsize',14);grid; set(gca,'fontsize',14);  %set(gca,'xlim',[0 delta_t*k],'ylim',[-0.05 0.05],'fontsize',14);
subplot(313); hold on; plot(Time,faults_sensor(3,1:k),'r','linewidth',2); plot(Time,p_af(1:k).*x_ob_filter_af(9,1:k),'b--','linewidth',2); ylabel('f_{\beta} (rad)','fontsize',14);grid; set(gca,'fontsize',14);  %set(gca,'xlim',[0 delta_t*k],'ylim',[-0.05 0.05],'fontsize',14);
h2=legend('True','DMAE'); set(h2,'color','white','edgecolor','black');
% h1=axes('position',[0.10 0.05 0.84 0.86],'fontsize',14);axis off;title(h1,'Probability-weighted Fault Estimation')
h1=axes('position',[0.50 0.0001 0.0001 0.0001],'fontsize',14); title('time (s)','fontsize',14)
% set(h2,'Box','off')


%-------------------  RMSE of fault estimation -------------------------
fprintf(1,'the RMSE of the state estimation is:\n%f\n%f\n%f\n ',...
    sqrt(mean((x_real(1,1:k)-x_ob_filter_af(1,1:k)).^2)), sqrt(mean((x_real(2,1:k)-x_ob_filter_af(2,1:k)).^2)),...
    sqrt(mean((x_real(3,1:k)-x_ob_filter_af(3,1:k)).^2)), sqrt(mean((x_real(4,1:k)-x_ob_filter_af(4,1:k)).^2)),...
    sqrt(mean((x_real(5,1:k)-x_ob_filter_af(5,1:k)).^2)), sqrt(mean((x_real(6,1:k)-x_ob_filter_af(6,1:k)).^2)));

% RMSE_f = sqrt(mean(error_f.^2'));
fprintf(1,'the RMSE of the fault estimation is:\n%f\n%f\n%f\n ',...
    sqrt(mean((faults_sensor(1,1:k) - p_af(1:k).*x_ob_filter_af(7,1:k)).^2)), ...
    sqrt(mean((faults_sensor(2,1:k) - p_af(1:k).*x_ob_filter_af(8,1:k)).^2)), ...
    sqrt(mean((faults_sensor(3,1:k) - p_af(1:k).*x_ob_filter_af(9,1:k)).^2)));
% fprintf(1,'the norm of root mean square error is:\n%f\n ',norm(RMSE_f(:)))
%-------------------  RMSE of fault estimation -------------------------


% state estimation
figure;
subplot(3,2,1);hold on; plot( Time,x_real(1,1:k),'r','linewidth',2); plot(Time,x_ob_filter_pdf2(1,1:k),'b--','linewidth',2);  ylabel('V (m/s)','fontsize',14);grid; set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);  %set(gcf,'color','none');set(gca,'color','none');% plot(P_sroot(1,1:k));plot( Time,-P_sroot(1,1:k));
% title('error of states, total')
subplot(3,2,2);hold on; plot( Time,x_real(2,1:k),'r','linewidth',2); plot(Time,x_ob_filter_pdf2(2,1:k),'b--','linewidth',2);  ylabel('\alpha (rad)','fontsize',14);grid;  set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14); % set(gcf,'color','none');set(gca,'color','none');% plot(P_sroot(2,1:k));plot( Time,-P_sroot(2,1:k));
subplot(3,2,3);hold on; plot( Time,x_real(3,1:k),'r','linewidth',2); plot(Time,x_ob_filter_pdf2(3,1:k),'b--','linewidth',2);  ylabel('\beta (rad)','fontsize',14);grid; set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);   %set(gcf,'color','none');set(gca,'color','none');%plot(P_sroot(3,1:k));plot( Time,-P_sroot(3,1:k));
subplot(3,2,4);hold on; plot( Time,x_real(4,1:k),'r','linewidth',2); plot(Time,x_ob_filter_pdf2(4,1:k),'b--','linewidth',2);  ylabel('\phi (rad)','fontsize',14);grid;  set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);  %set(gcf,'color','none');set(gca,'color','none');%plot(P_sroot(4,1:k));plot( Time,-P_sroot(4,1:k));
subplot(3,2,5);hold on; plot( Time,x_real(5,1:k),'r','linewidth',2); plot(Time,x_ob_filter_pdf2(5,1:k),'b--','linewidth',2);  ylabel('\theta (rad)','fontsize',14);grid; set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);  % set(gcf,'color','none');set(gca,'color','none');%plot(P_sroot(5,1:k));plot( Time,-P_sroot(5,1:k));
subplot(3,2,6);hold on; plot( Time,x_real(6,1:k),'r','linewidth',2); plot(Time,x_ob_filter_pdf2(6,1:k),'b--','linewidth',2);  ylabel('\psi (rad)','fontsize',14);grid;  set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);  %set(gcf,'color','none');set(gca,'color','none');%plot(P_sroot(6,1:k));plot( Time,-P_sroot(6,1:k));
% h=legend('ARafE','UKF'); set(h,'color','none','edgecolor','white');set(h,'Box','off')
h=legend('True','DMAE'); set(h,'Color','white','edgecolor','w'); % set(h,'Box','off')
h1=axes('position',[0.05 0.05 0.84 0.88],'fontsize',14);axis off;
% title(h1,'error of estimation')
h1=axes('position',[0.30 0.0001 0.0001 0.0001],'fontsize',14);%axis off;
title('time (s)','fontsize',14)
h1=axes('position',[0.80 0.0001 0.0001 0.0001],'fontsize',14);
title('time (s)','fontsize',14)


% state estimation error
figure;
subplot(3,2,1);hold on; plot( Time,x_real(1,1:k)-x_ob_filter_pdf2(1,1:k),'b','linewidth',2);  ylabel('\DeltaV (m/s)','fontsize',14);grid; set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);  set(gcf,'color','none');set(gca,'color','none');% plot(P_sroot(1,1:k));plot( Time,-P_sroot(1,1:k));
% title('error of states, total')
subplot(3,2,2);hold on; plot( Time,x_real(2,1:k)-x_ob_filter_pdf2(2,1:k),'b','linewidth',2);  ylabel('\Delta\alpha (rad)','fontsize',14);grid;  set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);  set(gcf,'color','none');set(gca,'color','none');% plot(P_sroot(2,1:k));plot( Time,-P_sroot(2,1:k));
subplot(3,2,3);hold on; plot( Time,x_real(3,1:k)-x_ob_filter_pdf2(3,1:k),'b','linewidth',2);  ylabel('\Delta\beta (rad)','fontsize',14);grid; set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);   set(gcf,'color','none');set(gca,'color','none');%plot(P_sroot(3,1:k));plot( Time,-P_sroot(3,1:k));
subplot(3,2,4);hold on; plot( Time,x_real(4,1:k)-x_ob_filter_pdf2(4,1:k),'b','linewidth',2);  ylabel('\Delta\phi (rad)','fontsize',14);grid;  set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);  set(gcf,'color','none');set(gca,'color','none');%plot(P_sroot(4,1:k));plot( Time,-P_sroot(4,1:k));
subplot(3,2,5);hold on; plot( Time,x_real(5,1:k)-x_ob_filter_pdf2(5,1:k),'b','linewidth',2);  ylabel('\Delta\theta (rad)','fontsize',14);grid; set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);   set(gcf,'color','none');set(gca,'color','none');%plot(P_sroot(5,1:k));plot( Time,-P_sroot(5,1:k));
subplot(3,2,6);hold on; plot( Time,x_real(6,1:k)-x_ob_filter_pdf2(6,1:k),'b','linewidth',2);  ylabel('\Delta\psi (rad)','fontsize',14);grid;  set(gca,'xlim',[0 delta_t*(k+1)], 'fontsize',14);  set(gcf,'color','none');set(gca,'color','none');%plot(P_sroot(6,1:k));plot( Time,-P_sroot(6,1:k));
h1=axes('position',[0.05 0.05 0.84 0.88],'fontsize',14);axis off;
% title(h1,'error of estimation')
h1=axes('position',[0.30 0.0001 0.0001 0.0001],'fontsize',14);%axis off;
title('time (s)','fontsize',14)
h1=axes('position',[0.80 0.0001 0.0001 0.0001],'fontsize',14);
title('time (s)','fontsize',14)

toc;