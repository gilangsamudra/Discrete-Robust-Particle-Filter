function total_error = maincost_func(i,F_hook,yini,setpoint)
% 3 steps
global con So y0 error1 error2 error3 waktu  

step1 = i;
So = F_hook(1:5);
if i == 1
    tspan = [0  waktu(step1)];
else 
    tspan = [0 waktu(step1)-waktu(step1-1)];
end
[t,y] = ode45('fun_ds_1dof_axial', tspan, yini);
y0    = y;
error1 = 1.0*abs(y0(end,2) - setpoint(step1));

step2 = i+1;
tspan = [0 waktu(step2)-waktu(step1)];
So = F_hook(6:10);
[t,y2] = ode45('fun_ds_1dof_axial', tspan, y0(end,:));
error2 = 1.0*abs(y2(end,2) - setpoint(step2));

step3 = i+2;
tspan = [0 waktu(step3)-waktu(step2)];
So = F_hook(11:15);
[t,y3] = ode45('fun_ds_1dof_axial', tspan, y2(end,:));
error3 = 1.0*abs(y3(end,2) - setpoint(step3));

total_error = error1 + error2 + error3;
con = [error1 error2 error3];
