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
function [inno,residual,error_x,Kk,Pk,Pzk,x_ob_filter,P_ob_filter,z_ob_filter]...
    =UKF_augmented_fault(mean_x0,P_x0,x_real,z_real,Qk,Rk,u,g,delta_t,namedas,w_c0,w_ci,w_m,dim_sys,k,Qd_est_M1)

%----------------start of augmented fault UKF----------------------
A_P=chol((dim_sys+namedas)*P_x0);% sqrt matrix 
A_P=A_P';
% sigma transform
mean_x0_copy=kron(mean_x0,ones(1,dim_sys));% copy the vector
sigma_x=[mean_x0  mean_x0_copy-A_P  mean_x0_copy+A_P];

V_g=mean_x0(1);
alpha_g=mean_x0(2);
beta_g=mean_x0(3);
phi_g=mean_x0(4);
theta_g=mean_x0(5);
psi_g=mean_x0(6);

% noise distribution matrix
G_noise(1,:)=[-cos(alpha_g)*cos(beta_g),-sin(beta_g),-sin(alpha_g)*cos(beta_g),...
    0,0,0];
G_noise(2,:)=[sin(alpha_g)/(V_g*cos(beta_g)),0,-cos(alpha_g)/(V_g*cos(beta_g)),...
    cos(alpha_g)*sin(beta_g)/cos(beta_g),-1,sin(alpha_g)*sin(beta_g)/cos(beta_g)];
G_noise(3,:)=[sin(beta_g)*cos(alpha_g)/V_g,-cos(beta_g)/V_g,sin(beta_g)*sin(alpha_g)/V_g,...
    -sin(alpha_g),0,cos(alpha_g)];
G_noise(4:6,:)=[    
    0  0  0     -1         -sin(phi_g)*tan(theta_g)        -cos(phi_g)*tan(theta_g);
    0  0  0      0              -cos(phi_g)                       sin(phi_g);
    0  0  0      0         -sin(phi_g)/cos(theta_g)         -cos(phi_g)/cos(theta_g);];

G_noise(7:9,:)=zeros(3,6); 


% time update
x_ob=model_sys_faulty(sigma_x,u(:,k),g(:,k),delta_t);
x_ob_mean=x_ob*w_m;

% variance
x_ob_mean_copy=kron(x_ob_mean,ones(1,2*dim_sys));

Qd_k=G_noise*Qk*G_noise';
% converted into the continuous noise covariance
Qd_k=Qd_k*delta_t;  

Pk=w_c0*(x_ob(:,1)-x_ob_mean)*(x_ob(:,1)-x_ob_mean)'...
    +w_ci*(x_ob(:,2:2*dim_sys+1)-x_ob_mean_copy)*(x_ob(:,2:2*dim_sys+1)-x_ob_mean_copy)'...
    +Qd_k;

%------------------------- use estimation of Qd -----------------------
Pk = Pk + Qd_est_M1;
%------------------------- use estimation of Qd -----------------------

Qfk = diag([zeros(1,9)]);
%
Pk = Pk + Qfk;


% nonlinear transform
z_ob=model_output_sensor_faults(x_ob);
% one step output
z_ob_mean=z_ob*w_m;

%%% measurement update
z_ob_mean_copy=kron(z_ob_mean,ones(1,2*dim_sys));
Pxz=w_c0*(x_ob(:,1)-x_ob_mean)*(z_ob(:,1)-z_ob_mean)'...
    +w_ci*(x_ob(:,2:2*dim_sys+1)-x_ob_mean_copy)*(z_ob(:,2:2*dim_sys+1)-z_ob_mean_copy)';

Pzk=Rk+w_c0*(z_ob(:,1)-z_ob_mean)*(z_ob(:,1)-z_ob_mean)'...
    +w_ci*(z_ob(:,2:2*dim_sys+1)-z_ob_mean_copy)*(z_ob(:,2:2*dim_sys+1)-z_ob_mean_copy)';

% Gain of Kalman Filter, update state and variance
Kk=Pxz/Pzk;

inno=z_real(:,k)-z_ob_mean;

% estimation 
x_ob_filter=x_ob_mean+Kk*inno;
P_ob_filter=Pk-Kk*Pzk*Kk';

% error & residual
error_x=x_real(:,k)-x_ob_filter(1:6,:);
z_ob_filter=model_output_sensor_faults(x_ob_filter);
residual=z_real(:,k)-z_ob_filter;




