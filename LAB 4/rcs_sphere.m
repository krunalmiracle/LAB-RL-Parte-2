function [] = rcs_sphere(i)
eps = 0.00001;
index = 0;
% kr limits are [0.05 - 15] ===> 300 points
for kr = 0.05:0.05:15
index = index + 1;
sphere_rcs = 0. + 0.*i;
f1 = 0. + 1.*i;
f2 = 1. + 0.*i;
m = 1.;
n = 0.;
q = -1.;
% initially set del to huge value
del =100000+100000*i;
while(abs(del) > eps)
q = -q;
n = n + 1;
m = m + 2;
del = (2.*n-1) * f2 / kr-f1;
f1 = f2;
f2 = del;
del = q * m /(f2 * (kr * f1 - n * f2));
sphere_rcs = sphere_rcs + del;
end
rcs(index) = abs(sphere_rcs);
sphere_rcsdb(index) = 10. * log10(rcs(index));
end
figure; set(gcf,'Color','w');
n=0.05:.05:15;
plot (n,rcs,'k');
set (gca,'xtick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]);
xlabel ('Sphere circumference in wavelengths (2{\pi}r/{\lambda})');
ylabel ('Normalized sphere RCS [adim.]');
grid;
figure; set(gcf,'Color','w');
plot (n,sphere_rcsdb,'k');
set (gca,'xtick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]);
xlabel ('Sphere circumference in wavelengths (2{\pi}r/{\lambda})');
ylabel ('Normalized sphere RCS [dB]');
grid;
figure; set(gcf,'Color','w');
semilogx (n,sphere_rcsdb,'k'); grid on;
xlabel ('Sphere circumference in wavelengths (2{\pi}r/{\lambda})');
ylabel ('Normalized sphere RCS [dB]');
end
