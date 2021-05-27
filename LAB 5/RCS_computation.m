function [monostatic_RCS]=RCS_computation(RCS_scatters,scatters_rect_coord,radar_sphere_coord,carrier_freq)
    %coord esfèriques dels radars
    r=radar_sphere_coord(:,1);
    theta=radar_sphere_coord(:,2);
    phi=radar_sphere_coord(:,3);
    %coord cartesianes dels radars
    x_radar=r.*sin(theta).*cos(phi);
    y_radar=r.*sin(theta).*sin(phi);
    z_radar=r.*cos(theta);
    %coord cartesianes dels scatters
    x_scatter=scatters_rect_coord(:,1);
    y_scatter=scatters_rect_coord(:,2);
    z_scatter=scatters_rect_coord(:,3);
    c=3*10^8;
    lambda=c/carrier_freq;
    RCS_scatters=10.^(RCS_scatters/10);
    %per cada radar-->tots els scatters
    for k=1:length(x_radar)
    for j=1:length(x_scatter)
    r_n(k,j)=sqrt((x_scatter(j)-x_radar(k)).^2+(y_scatter(j)-
    y_radar(k)).^2+(z_scatter(j)-z_radar(k)).^2);
    phi_n(k,j)=2*2*pi/lambda*r_n(k,j);
    end
    monostatic_RCS(k)=abs(sum(sqrt(RCS_scatters).*exp(1i*phi_n(k,:)))).^2;
    end
end
