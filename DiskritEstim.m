clear;clc;
tic
load('Input Accel.mat');
load('Viscous Force Approximation.mat');
Hook = [0 Trajectory(2:end)];

%% Parameter Particle Filter
process = 0.00125;
measure = 0.1;
w = process*ones(2,10);
z = measure*ones(1,1);
num_members=100;
alp_k = 0.1;
eps   = 0.5;
p1 = 10;
fre = 100;
LER = 0.05;
stdnom = 0.25;
%% DRAG FORCES
Loaded_Data = readtable('Drag Calculation.xlsx'); % in table format
drag        = Loaded_Data{:,5};                      % drag forces (lbf)
drag = drag(~isnan(drag));

% Panjang in meter
Total = 5500;
ns    = 4;
LengthS = Total/ns;
LengthSNom = Total/ns - (LER*LengthS);
mod_el = 206842718795.3;          % Elastic modulus, N/m^2

% String Table, inches
ODm = 5 * 0.0254; %in meter now
ODj = 6.63 * 0.0254;
ODc = 5 * 0.0254;
ODb = 6.85 * 0.0254;
IDm = 4.28 * 0.0254;
IDj = 3.25 * 0.0254;
IDc = 3.19 * 0.0254;
IDb = 3.89 * 0.0254;

% in SI unit
pipe_den = 7800.001722079;        %Pipe density, in kg/m^3

% Drillstring Parameters
Stiff_bha = mod_el*((ODb^2 - IDb^2)*pi/4)/(0.088*LengthS);
Stiff_coll= mod_el*((ODc^2 - IDc^2)*pi/4)/(0.911*LengthS);
Stiff_main= mod_el*((ODm^2 - IDm^2)*pi/4)/(0.95*LengthS);
Stiff_sec = mod_el*((ODj^2 - IDj^2)*pi/4)/(0.05*LengthS);

stiffness1 = Stiff_main*Stiff_sec/(Stiff_main+Stiff_sec);     %stiffness of drill pipe
stiffness2 = stiffness1;
stiffness3 = stiffness1;
stiffness4 = Stiff_bha*Stiff_coll/(Stiff_bha+Stiff_coll);
stiffness = [stiffness1 stiffness2 stiffness3 stiffness4];

mass_top  = 20000;
mass_bha  = ((ODb^2 - IDb^2)*pi/4)*(0.088*LengthS)*pipe_den;        
mass_coll = ((ODc^2 - IDc^2)*pi/4)*(0.911*LengthS)*pipe_den;
mass_main = ((ODm^2 - IDm^2)*pi/4)*(0.95*LengthS)*pipe_den;
mass_sec  = ((ODj^2 - IDj^2)*pi/4)*(0.05*LengthS)*pipe_den;

mass1 = mass_main+mass_sec; %mass of drill pipe
mass2 = mass1;
mass3 = mass1;
mass4 = mass_bha+mass_coll;
mass  = [mass_top mass1 mass2 mass3 mass4];


r1 = stiffness1/mass_top;
r12= stiffness1/mass1;
r2 = stiffness2/mass1;
r23= stiffness2/mass2;
r3 = stiffness3/mass2;
r34= stiffness3/mass3;
r4 = stiffness4/mass3;
r45= stiffness4/mass4;
r = [r1 r12 r2 r23 r3 r34 r4 r45];

% True Plant
A = [0 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 0 1 0 0 0 0 0 0
    r12 0 (-r12-r2) -pp1(1)/mass1 r2 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0
    0 0 r23 0 (-r23-r3) -pp2(1)/mass2 r3 0 0 0
    0 0 0 0 0 0 0 1 0 0
    0 0 0 0 r34 0 (-r34-r4) -pp3(1)/mass3 r4 0
    0 0 0 0 0 0 0 0 0 1
    0 0 0 0 0 0 r45 0 -r45 -pp4(1)/mass4];

B = [0 0 0 0 0 0 0 0 0 0
    1 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 -1/mass1 0 0 0 -1/mass1 0 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 -1/mass2 0 0 0 -1/mass2 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 0 -1/mass3 0 0 0 -1/mass3 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 0 0 -1/mass4 0 0 0 -1/mass4 0];
 
C = eye(10);
D = zeros(10,10);

%% Drillstring Parameters Nominal
Stiff_bhaNom = mod_el*((ODb^2 - IDb^2)*pi/4)/(0.088*LengthSNom);
Stiff_collNom= mod_el*((ODc^2 - IDc^2)*pi/4)/(0.911*LengthSNom);
Stiff_mainNom= mod_el*((ODm^2 - IDm^2)*pi/4)/(0.95*LengthSNom);
Stiff_secNom = mod_el*((ODj^2 - IDj^2)*pi/4)/(0.05*LengthSNom);

stiffness1Nom = Stiff_mainNom*Stiff_secNom/(Stiff_mainNom+Stiff_secNom);     %stiffness of drill pipe
stiffness2Nom = stiffness1Nom;
stiffness3Nom = stiffness1Nom;
stiffness4Nom = Stiff_bhaNom*Stiff_collNom/(Stiff_bhaNom+Stiff_collNom);
stiffnessNom = [stiffness1Nom stiffness2Nom stiffness3Nom stiffness4Nom];

mass_bhaNom  = ((ODb^2 - IDb^2)*pi/4)*(0.088*LengthSNom)*pipe_den;        
mass_collNom = ((ODc^2 - IDc^2)*pi/4)*(0.911*LengthSNom)*pipe_den;
mass_mainNom = ((ODm^2 - IDm^2)*pi/4)*(0.95*LengthSNom)*pipe_den;
mass_secNom  = ((ODj^2 - IDj^2)*pi/4)*(0.05*LengthSNom)*pipe_den;

mass1Nom = mass_mainNom+mass_secNom; %mass of drill pipe
mass2Nom = mass1Nom;
mass3Nom = mass1Nom;
mass4Nom = mass_bhaNom+mass_collNom;
massNom  = [mass_top mass1Nom mass2Nom mass3Nom mass4Nom];

r1nom = stiffness1Nom/mass_top;
r12nom= stiffness1Nom/mass1Nom;
r2nom = stiffness2Nom/mass1Nom;
r23nom= stiffness2Nom/mass2Nom;
r3nom = stiffness3Nom/mass2Nom;
r34nom= stiffness3Nom/mass3Nom;
r4nom = stiffness4Nom/mass3Nom;
r45nom= stiffness4Nom/mass4Nom;
rNom = [r1nom r12nom r2nom r23nom r3nom r34nom r4nom r45nom];

% Modeled Plant (nominal)
Anom = [0 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 0 1 0 0 0 0 0 0
    r12nom 0 (-r12nom-r2nom) -pp1(1)/mass1 r2nom 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0
    0 0 r23nom 0 (-r23nom-r3nom) -pp2(1)/mass2 r3nom 0 0 0
    0 0 0 0 0 0 0 1 0 0
    0 0 0 0 r34nom 0 (-r34nom-r4nom) -pp3(1)/mass3 r4nom 0
    0 0 0 0 0 0 0 0 0 1
    0 0 0 0 0 0 r45nom 0 -r45nom -pp4(1)/mass4];

% Uncertainty Matrix
A_bar = [0 0 0 0 0 0 0 0 0 0
         0 0 0 0 0 0 0 0 0 0
         0 0 0 0 0 0 0 0 0 0
         r12nom*alp_k 0 (-r12nom-r2nom)*alp_k 0 r2nom*alp_k 0 0 0 0 0
         0 0 0 0 0 0 0 0 0 0
         0 0 r23nom*alp_k 0 (-r23nom-r3nom)*alp_k 0 r3nom*alp_k 0 0 0
         0 0 0 0 0 0 0 0 0 0
         0 0 0 0 r34nom*alp_k 0 (-r34nom-r4nom)*alp_k 0 r4nom*alp_k 0
         0 0 0 0 0 0 0 0 0 0
         0 0 0 0 0 0 r45nom*alp_k 0 -r45nom*alp_k 0];

% untuk simulasi
t = 0:(1/fre):30;
drag1 = drag(1)*ones(numel(t),1);
drag2 = drag(2)*ones(numel(t),1);
drag3 = drag(3)*ones(numel(t),1);
drag4 = drag(4)*ones(numel(t),1);
vis1  = pp1(2)*ones(numel(t),1);
vis2  = pp2(2)*ones(numel(t),1);
vis3  = pp3(2)*ones(numel(t),1);
vis4  = pp4(2)*ones(numel(t),1);

hook2 = [Hook(1)];
for i=1:numel(Hook)-1
    interpol = linspace(Hook(i),Hook(i+1),fre);
    hook2 = [hook2 interpol];
end
input1 = [hook2' drag1 drag2 drag3 drag4 vis1 vis2 vis3 vis4 zeros(numel(t),1)];
input2 = [(hook2'-(0.0*hook2')) drag1 drag2 drag3 drag4 vis1 vis2 vis3 vis4 zeros(numel(t),1)];


contrue = ss(A,B,C,D);
connomi = ss(Anom,B,C,D);
disktrue = c2d(contrue,(1/fre));
disknomi = c2d(connomi,(1/fre));
xo = zeros(10,1);

% PF Routine
[xtrue,xestm]=diskrit_PF(A_bar,input1,input2,t,xo,disktrue,disknomi,p1,w,z,num_members,eps,stdnom);
error = abs(xtrue-xestm);
rataerr = mean(error(:,10));

% Plotting
figure()
plot(t(1:end-1),xtrue(:,10),'LineWidth',1);
hold on;
plot(t(1:end-1),xestm(:,10),'LineWidth',1);
set(gca,'FontSize',15, 'xlim',[t(1) t(end)+1],'ylim',[min(xtrue(:,10))-1 max(xtrue(:,10))+1])
xlabel('Time (s)','FontSize', 24)
ylabel('BHA Velocity (m/s)','FontSize', 24)
legend('Real Velocity','Estimated Velocity')
text(5,-2,['Number of Ensemble = ',num2str(num_members)],'fontsize',13);
text(5,-2.5,['Uncertainty = ',num2str(alp_k*100),'%'],'fontsize',13);
text(5,-3,['Frequency Sampling = ',num2str(fre),'Hz'],'fontsize',13);
text(5,-3.5,['Length Error = ',num2str(LER*100),'%'],'fontsize',13);
text(5,-4,['Epsilon = ',num2str(eps)],'fontsize',13);
text(5,-4.5,['Normal Std.Dev. = ',num2str(stdnom*100),'%'],'fontsize',13);
text(5,-5,['Process Noise Std.Dev. = ',num2str(process)],'fontsize',13);
text(5,-5.5,['Measurement Noise Std.Dev. = ',num2str(measure)],'fontsize',13);
hold off
% 
% figure()
% plot(t(1:end-1),error(:,10),'LineWidth',1);
% set(gca,'FontSize',24, 'xlim',[t(1) t(end)+1],'ylim',[min(xtrue(:,10))-1 max(xtrue(:,10))+1])
% xlabel('Time (s)','FontSize', 24)
% ylabel('Error','FontSize', 24)
% text(5,-3,['Number of Ensemble = ',num2str(num_members)],'fontsize',18);
% text(5,-4,['Uncertainty = ',num2str(alp_k*100),'%'],'fontsize',18);
% text(5,-5,'Nominal Model only','fontsize',18);
% 
% figure()
% plot(t(1:end-1),xtrue(:,2),'LineWidth',1);
% hold on;
% plot(t(1:end-1),xestm(:,2),'LineWidth',1);
% set(gca,'FontSize',24, 'xlim',[t(1) t(end)+1],'ylim',[min(xestm(:,2))-1 max(xestm(:,2))+1])
% xlabel('Time (s)','FontSize', 24)
% ylabel('BHA Velocity (m/s)','FontSize', 24)
% legend('Real Velocity','Estimated Velocity')
% text(5,-3,['Number of Ensemble = ',num2str(num_members)],'fontsize',18);
% text(5,-4,['Uncertainty = ',num2str(alp_k*100),'%'],'fontsize',18);
% text(5,-5,'Nominal Model only','fontsize',18);
% hold off
toc