% subfunctions 
% augmented states of the observer
function x_ob=model_sys_faulty(x,u,g,delta_t)

Vtas=x(1,:);alpha=x(2,:);beta=x(3,:);
phi=x(4,:);theta=x(5,:);psi=x(6,:);

% ignore these, these are non-zero if you consider IMU FDD
la_x=0;la_y=0;la_z=0;
la_p=0;la_q=0;la_r=0;

Axm=u(1,:);Aym=u(2,:);Azm=u(3,:);
pm=u(4,:);qm=u(5,:);rm=u(6,:);   
    
% the fault of the measurements
delta_V=x(7,:);
delta_a=x(8,:);
delta_b=x(9,:);


% note that this function allows multiple-dimension x
[~,den]=size(x);

x_ob=[  Vtas+delta_t.*(cos(alpha).*cos(beta).*((Axm-la_x)-g.*sin(theta)+(rm-la_r).*Vtas.*sin(beta)-(qm-la_q).*Vtas.*sin(alpha).*cos(beta))...
    +sin(beta).*((Aym-la_y)+g.*cos(theta).*sin(phi)+(pm-la_p).*Vtas.*sin(alpha).*cos(beta)-(rm-la_r).*Vtas.*cos(alpha).*cos(beta))...
    +sin(alpha).*cos(beta).*((Azm-la_z)+g.*cos(theta).*cos(phi)+(qm-la_q).*Vtas.*cos(alpha).*cos(beta)-(pm-la_p).*Vtas.*sin(beta)));
        alpha+delta_t.*(cos(alpha)./(Vtas.*cos(beta)).*((Azm-la_z)+g.*cos(theta).*cos(phi)+(qm-la_q).*Vtas.*cos(alpha).*cos(beta)-(pm-la_p).*Vtas.*sin(beta))...
    -sin(alpha)./(Vtas.*cos(beta)).*((Axm-la_x)-g.*sin(theta)+(rm-la_r).*Vtas.*sin(beta)-(qm-la_q).*Vtas.*sin(alpha).*cos(beta)));
        beta+delta_t.*(1./(Vtas.*cos(beta)).*((Aym-la_y)+g.*cos(theta).*sin(phi)+(pm-la_p).*Vtas.*sin(alpha).*cos(beta)-(rm-la_r).*Vtas.*cos(alpha).*cos(beta))...
    -sin(beta)./(Vtas.*cos(beta)).*(cos(alpha).*cos(beta).*((Axm-la_x)-g.*sin(theta)+(rm-la_r).*Vtas.*sin(beta)-(qm-la_q).*Vtas.*sin(alpha).*cos(beta))...
    +sin(beta).*((Aym-la_y)+g.*cos(theta).*sin(phi)+(pm-la_p).*Vtas.*sin(alpha).*cos(beta)-(rm-la_r).*Vtas.*cos(alpha).*cos(beta))...
    +sin(alpha).*cos(beta).*((Azm-la_z)+g.*cos(theta).*cos(phi)+(qm-la_q).*Vtas.*cos(alpha).*cos(beta)-(pm-la_p).*Vtas.*sin(beta))));
        phi+delta_t.*((pm-la_p)+(qm-la_q).*sin(phi).*tan(theta)+(rm-la_r).*cos(phi).*tan(theta));
        theta+delta_t.*((qm-la_q).*cos(phi)-(rm-la_r).*sin(phi));
        psi+delta_t.*((qm-la_q).*sin(phi)./cos(theta)+(rm-la_r).*cos(phi)./cos(theta));
        delta_V;
        delta_a;
        delta_b;];

