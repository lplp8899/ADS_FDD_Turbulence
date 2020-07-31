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
function G_noise=G_Ja(x)

V_g=x(1,:);alpha_g=x(2,:);beta_g=x(3,:);
phi_g=x(4,:);theta_g=x(5,:);psi_g=x(6,:);

%
G_noise(1,:)=[-cos(alpha_g)*cos(beta_g),-sin(beta_g),-sin(alpha_g)*cos(beta_g),...
    0,0,0];
G_noise(2,:)=[sin(alpha_g)/(V_g*cos(beta_g)),0,-cos(alpha_g)/(V_g*cos(beta_g)),...
    cos(alpha_g)*tan(beta_g),-1,sin(alpha_g)*tan(beta_g)];
G_noise(3,:)=[sin(beta_g)*cos(alpha_g)/V_g,-cos(beta_g)/V_g,sin(beta_g)*sin(alpha_g)/V_g,...
    -sin(alpha_g),0,cos(alpha_g)];
G_noise(4:6,:)=[    
    0  0  0     -1         -sin(phi_g)*tan(theta_g)        -cos(phi_g)*tan(theta_g);
    0  0  0      0              -cos(phi_g)                       sin(phi_g);
    0  0  0      0         -sin(phi_g)/cos(theta_g)         -cos(phi_g)/cos(theta_g);];



