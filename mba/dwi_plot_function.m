function []= dwi_plot_function(dwi_vector, x_dwi, y_dwi, z_dwi)
    % surf spherical functions
    the_cnt = 48;
    phi_cnt = 96;
    theta3 = linspace(0*pi/180,180*pi/180,the_cnt);
    phi3 = linspace(0*pi/180,360*pi/180,phi_cnt);
    [theta, phi] = meshgrid(theta3,phi3);

    F = scatteredInterpolant(x_dwi, y_dwi, z_dwi, dwi_vector);
    F.Method = 'natural';

    x = 1.*sin(theta).*cos(phi);
    y = 1.*sin(theta).*sin(phi);
    z = 1.*cos(theta);
    t1 = reshape(x, [96*48 1]);
    t2 = reshape(y, [96*48 1]);
    t3 = reshape(z, [96*48 1]);

    dwi_interpolated = F(t1, t2, t3);
    T = reshape(dwi_interpolated,[96 48]);
    
    x2 = T.*sin(theta).*cos(phi);
    y2 = T.*sin(theta).*sin(phi);
    z2 = T.*cos(theta);

    surf(x2,y2,z2,T);
    shading interp
    axis off
end
