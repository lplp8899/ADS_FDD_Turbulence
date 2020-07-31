function y_ob=model_output_sensor_faults(x)

delta_V=x(7,:);
delta_a=x(8,:);
delta_b=x(9,:);

% output model with sensor faults
   
y_ob=[ x(1,:)+delta_V;
       x(2,:)+delta_a;
       x(3,:)+delta_b;
       x(4,:);
       x(5,:);
       x(6,:);];   
   
