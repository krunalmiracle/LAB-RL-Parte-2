clear all;
clear all figures;
%% 5.1 RCS of an sphere --> dependency: frequency of the radar
%% 5.2 RCS of a circular flat plate --> dependency: aspect angle
%INPUTS: frequency [Hz], radius [m]
%OUTPUT: aspect angle [º], RCS [dB]
r=10000;
freq=10;
[vaspect_deg, rcs_dB] = rcs_circ_plate (r, freq);
lambda=3e8/freq;
area=pi*r^2;
%comprovem que per theta=0º-->RCS=(4pi*Area^2)/(lambda^2)
normal_incidence_RCS=10*log10((4*pi*area^2)/(lambda^2));
%% 5.4 Two spheres RCS
%==============RCS COMPUTATION FUNCTION====================
%-----I N P U T S------
%RCS_scatters = Column vector of RCSs of spheres [dBsm=dBm2]
%scatters_rect_coord = Matrix having the rectangular coordinates of the scatterers positions:
% x-coordinates in first column [m]
% y-coordinates in second column [m]
% z-coordinates in third column [m]
%radar_sphere_coord = Matrix having the spherical coordinates of the different radar locations
% r-coordinates in first column [m]
% theta-coordinates in second column [rad]
% phi-coordinates in third column [rad]
%carrier_freq=Carrier frequency [Hz]
%-----O U T P U T S------
% monostatic_RCS = Column vector of monostatic RCSs seen from the different radar locations [m^2]
%==============FOR ONE VALUE OF L==========================
%--> i n p u t s a j u s t a b l e s
L=sqrt(3)/2; %-->electrical separation between scatters
carrier_freq=1.5*10^9;
RCS_scatters=[0,0, 0, 0, 16]; %0dBsm-->1m^2 isotropical scatters
%coordinates of scatters // x y z
scatters_rect_coord=[0 0 0;sqrt(3)*L 0 0;0 L 0;L sqrt(3)*L 0;sqrt(3)*L/2 L/2 0];
%coordinates of radars
r_gran=10^6;
if(r_gran<(2*L^2)*carrier_freq/(3e8)) % condition to fulfill
    r_gran<(2*L^2/lambda)
    disp('Distance between scatters and radar stations too small.');
end
theta=pi/2; %-->horizontal plane
phi=0:pi/1000:2*pi; %--> all angles of the horizontal plane
%--> e x e c u c i ó
for i=1:length(phi)
    radar_sphere_coord(i,:)=[r_gran theta phi(i)];
end
monostatic_RCS=RCS_computation(RCS_scatters,scatters_rect_coord,radar_sphere_coord,carrier_freq);
figure(1);
hold on;
plot(phi,10*log10(monostatic_RCS)); %RCS [dBsm]
hold on;
xlabel('Phi angle (rad)');
hold on;
ylabel('RCS (dBsm)');
hold on;
str = sprintf('Resulting RCS for five spheres separated %0.2fm',L);
hold on;
title(str);
hold off;
figure(2);
polarscatter(phi,monostatic_RCS); %RCS [m^2]
str = sprintf('Resulting RCS for five spheres (1m^2) separated %0.2fm',L);
title(str);
%===============FOR DIFFERENT VALUES OF L=====================
%--> i n p u t s a j u s t a b l e s
% L=[1 0.1 0.01];
j=1;
carrier_freq=1.5*10^9;
while(j<=length(L))
%     RCS_scatters=[0,0]; %0dBsm-->1m^2
%     scatters_rect_coord=[-L(j)/2 0 0; L(j)/2 0 0];
    r_gran=10^6;
    theta=pi/2;
    phi=0:pi/1000:2*pi;
    %--> e x e c u c i ó
    for i=1:length(phi)
        radar_sphere_coord(i,:)=[r_gran theta phi(i)];
    end
    monostatic_RCS=RCS_computation(RCS_scatters,scatters_rect_coord,radar_sphere_coord,carrier_freq);
    figure(3);
    subplot(length(L),2,j+(j-1))
    plot(phi,10*log10(monostatic_RCS)); %RCS [dBsm]
    xlabel('Phi angle (rad)');
    ylabel('RCS (dBsm)');
    subplot(length(L),2,(j+1)+(j-1))
    polarplot(phi,monostatic_RCS); %RCS [m^2]
    str = sprintf('Resulting RCS for four spheres with L separated %0.2f m',L(j));
    title(str);
    j=j+1;
end
%% 5.5 Swerling model statistics I and II
%==============RCS COMPUTATION FUNCTION====================
%-----I N P U T S------
%RCS_scatters = Column vector of RCSs of spheres [dBsm=dBm2]
%scatters_rect_coord = Matrix having the rectangular coordinates of the scatterers positions:
% x-coordinates in first column [m]
% y-coordinates in second column [m]
% z-coordinates in third column [m]
%radar_sphere_coord = Matrix having the spherical coordinates of the different radar locations
% r-coordinates in first column [m]
% theta-coordinates in second column [rad]
% phi-coordinates in third column [rad]
%carrier_freq=Carrier frequency [Hz]
%-----O U T P U T S------
% monostatic_RCS = Column vector of monostatic RCSs seen from the different radar locations [m^2]
%--> i n p u t s a j u s t a b l e s
% L=1;
carrier_freq=3*10^9;
% num_scatters=5; %10 scatters in the scenario
% RCS_scatters=ones(1,num_scatters)*0.0000001; %dBsm RCSs close to 0 dBsm;
%coordinates of scatters DISTRIBUTED RANDOMLY
% circle_radius=10; %circle of 10m
% r_scatter=0+circle_radius*rand(1,num_scatters);
% theta_scatter=ones(1,num_scatters)*pi/2;
% phi_scatter=0+2*pi*rand(1,num_scatters);
% x_scatter=r_scatter.*sin(theta_scatter).*cos(phi_scatter);
% y_scatter=r_scatter.*sin(theta_scatter).*sin(phi_scatter);
% z_scatter=r_scatter.*cos(theta_scatter);
% scatters_rect_coord=[x_scatter; y_scatter; z_scatter]';
%coordinates of radar positions
r_gran=1000; %1km
theta=pi/2; %-->horizontal plane
phi=0:pi/1000:2*pi;%--> all angles of the horizontal plane
%--> e x e c u c i ó
for i=1:length(phi)
    radar_sphere_coord(i,:)=[r_gran theta phi(i)];
end
monostatic_RCS=RCS_computation(RCS_scatters,scatters_rect_coord,radar_sphere_coord,carrier_freq);
figure(1);
histogram(monostatic_RCS,'Normalization','pdf'); %RCS [m^2]
data_hist_1=histcounts(monostatic_RCS,'Normalization','pdf');
ylabel('RCS probability density function');
xlabel('RCS (m^2)');
title(sprintf('PDF of RCS resulting from %g spheres with 0dBsm and 16dbsm distributed ',5));
figure(2);
subplot(2,1,1)
plot(phi,10*log10(monostatic_RCS)); %RCS [dBsm]
xlabel('Phi angle (rad)');
ylabel('RCS (dBsm)');
subplot(2,1,2)
polarplot(phi,monostatic_RCS); %RCS [m^2]
subtitle(sprintf('Resulting RCS for %g spheres with 0dBsm and 16dbsm distributed ',5));
%% 5.6 Swerling model statistics III and IV
%==============RCS COMPUTATION FUNCTION====================
%-----I N P U T S------
%RCS_scatters = Column vector of RCSs of spheres [dBsm=dBm2]
%scatters_rect_coord = Matrix having the rectangular coordinates of the scatterers positions:
% x-coordinates in first column [m]
% y-coordinates in second column [m]
% z-coordinates in third column [m]
%radar_sphere_coord = Matrix having the spherical coordinates of the different radar locations
% r-coordinates in first column [m]
% theta-coordinates in second column [rad]
% phi-coordinates in third column [rad]
%carrier_freq=Carrier frequency [Hz]
%-----O U T P U T S------
% monostatic_RCS = Column vector of monostatic RCSs seen from the different radar locations [m^2]
%--> i n p u t s a j u s t a b l e s
L=1;
carrier_freq=3*10^9;
% num_scatters=5; %10 scatters in the scenario
% RCS_scatters=ones(1,num_scatters)*0.0000001; %dBsm RCSs close to 0 dBsm;
% RCS_scatters(8)=RCS_scatters(8)+17; %one scatter 17dB greater
%coordinates of scatters DISTRIBUTED RANDOMLY
% circle_radius=10; %circle of 10m
% r_scatter=0+circle_radius*rand(1,num_scatters);
% theta_scatter=ones(1,num_scatters)*pi/2;
% phi_scatter=0+2*pi*rand(1,num_scatters);
% x_scatter=r_scatter.*sin(theta_scatter).*cos(phi_scatter);
% y_scatter=r_scatter.*sin(theta_scatter).*sin(phi_scatter);
% z_scatter=r_scatter.*cos(theta_scatter);
% scatters_rect_coord=[x_scatter; y_scatter; z_scatter]';
%coordinates of radar positions
r_gran=1000; %1km
theta=pi/2; %-->horizontal plane
phi=0:pi/1000:2*pi;%--> all angles of the horizontal plane
%--> e x e c u c i ó
for i=1:length(phi)
    radar_sphere_coord(i,:)=[r_gran theta phi(i)];
end
monostatic_RCS=RCS_computation(RCS_scatters,scatters_rect_coord,radar_sphere_coord,carrier_freq);
figure(3);
histogram(monostatic_RCS,'Normalization','pdf'); %RCS [m^2]
data_hist_1=histcounts(monostatic_RCS,'Normalization','pdf');
ylabel('RCS probability density function');
xlabel('RCS (m^2)');
title(sprintf('PDF of RCS resulting from %g spheres with 0dBsm except 1 with 16 distributed ',5));
figure(4);
subplot(2,1,1)
plot(phi,10*log10(monostatic_RCS)); %RCS [dBsm]
xlabel('Phi angle (rad)');
ylabel('RCS (dBsm)');
subplot(2,1,2)
polarplot(phi,monostatic_RCS); %RCS [m^2]
subtitle(sprintf('Resulting RCS for %g spheres with 0dBsm except 1 with 16 distributed',5));