clear;clc;%close all
tic
global m c k So y0 waktu input  
%% Speed Trajectory
fre = 100;
setpoint = [linspace(0,0.51,fre) ones(1,fre)*0.51 linspace(0.51,0,fre)];
waktu = linspace(0.5,30,numel(setpoint));
%% Dapeng's Case
iCase = 1;   % select case for mass (1 to 5)

% pipe properties
rho_pipe =   30;              % kg/m³, mass/length of pipe
rho_bha  =  200;              % kg/m³, mass/length of bha
L_pipe   = 5000;              % m, length of pipe
L_bha    =  100;              % m, length of bha
m_pipe   = rho_pipe * L_pipe; % kg, mass of pipe
m_bha    = rho_bha * L_bha;   % kg, mass of bha
m_bit    =  200;              % kg, mass of bit
D_pipe   = 5.0 * 0.0254;      % inch --> m, outer diameter of pipe
D_hole   = 8.5 * 0.0254;      % inch --> m, diameter of bore hole
H        = 1/2*(D_hole-D_pipe); % m, average annular clearance
k        = 115000;            % N/m; stiffness of drill string
c        =  30000;            % N·s/m; damping constant

% select case for mass
switch iCase
    case 1
        m    = (m_bha + m_pipe + m_bit)*1; % kg;
    case 2
        m    = m_bha + 1/3 * m_pipe; % kg;
    case 3
        m    = m_bha + m_bit; % kg;
    case 4
        m    = m_bha; % kg;
    case 5
        m    = m_bit; % kg;
end

%% Optimizer
% Initial F hook
x0 = ones(1,5*3); % Initial accelerations untuk t0, t0+1, t0+2

% Lower & Upper Bounds
lb = ones(1,5*3)*0;
ub = ones(1,5*3)*8e6;

options = optimoptions(@fmincon,...
    'MaxIter', 500, ...
    'MaxFunEvals', 10000,...
    'TolFun',1e-6,...
    'Display', 'iter', ...
    'Algorithm', 'active-set');

input   = [];
keluaran= [];
kecepatan = [];
yini = [0;0];
for i = 1:numel(setpoint)-2
    disp(['step ke ',num2str(i)])
    [F_hook,evaluasi] = fmincon(@(F_hook)maincost_func(i,F_hook,yini,setpoint),...
        x0,[],[],[],[],lb,ub,'nonlicon2',options);
    input = [input F_hook(1:5)];
    keluaran = [keluaran evaluasi(1)];
    
    step1 = i;
    So = F_hook(1:5);
    if i == 1
        tspan = [0  waktu(step1)];
    else
        tspan = [0 waktu(step1)-waktu(step1-1)];
    end
    [~,y] = ode45('fun_ds_1dof_axial', tspan, y0(end,:));
    yini = y(end,:);
    kecepatan= [kecepatan y(end,2)];
end
selisih = setpoint(end-2) - kecepatan;

figure()
plot(waktu,setpoint,'--','linewidth',3)
hold on
plot(waktu(1:end-2), kecepatan,'linewidth',3)
hold off
title('Speed Trajectory')
xlabel('Time')
ylabel('Speed')
legend('Speed Setpoint','MPC')
set(gca,'fontsize',15)
toc