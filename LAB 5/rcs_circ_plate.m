function [vaspect_deg, rcs_dB] = rcs_circ_plate (r, freq)

% This function calculates and plots the RCS of a circular flat plate of radius r.
% (c) B.H. Mahafza, Radar Systems Analysis and Design Using MATLAB,
%       Chapman & Hall/CRC, 2005
% Modified by J.M. González-Arbesú

eps = 0.000001;
cdeg2rad=pi/180;

% Compute wavelength
lambda=3e8/freq;
index=0;
vaspect_deg=[0:.1:180];
vaspect_rad=cdeg2rad*vaspect_deg;
vx=(4*pi*r/lambda)*sin(vaspect_rad);
vval1=4*pi^3*r^4/lambda^2;
vval2=2*besselj(1,vx)./vx;
rcs_po2=vval1*(vval2.*cos(vaspect_rad)).^2+eps;
% vval1m=lambda*r;
% vval2m=8*pi*sin(vaspect_rad).*(tan(vaspect_rad).^2);
% rcs_mu2=vval1m./vval2m+eps;
[vindex]=find(or(vaspect_deg==0,vaspect_deg==180));
rcs_po2(vindex)=(4*pi^3*r^4/lambda^2)+eps;
% rcs_mu2(vindex)=rcs_po2(vindex);

rcsdb_po2=10*log10(rcs_po2);
% rcsdb_mu2=10*log10(rcs_mu2);

figure; set(gcf,'Color','w');
plot(vaspect_deg,rcsdb_po2,'b'); hold on; grid on;
% plot(vaspect_deg,rcsdb_mu2,'b--');
xlabel ('Aspect angle [º]');
ylabel ('RCS [dBsm]');
rcs_dB=[rcsdb_po2];
