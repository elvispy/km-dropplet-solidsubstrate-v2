function rmax = maximum_contact_radius(amplitudes)
    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
    %collectdnPl = PROBLEM_CONSTANTS.collectdnPl;
    %collectd2nPl = PROBLEM_CONSTANTS.collectd2nPl;
    order = length(amplitudes);
    zeta = zeta_generator(amplitudes); 
    drdtheta = @(theta) cos(theta) * (1 + zeta(theta)) - ...
        sin(theta)^2 * sum(times(amplitudes, collectdnPl(order, cos(theta))));

    dr2dtheta2 = @(theta) - sin(theta) * (1 + zeta(theta)) - 2 * cos(theta) * sin(theta) * ...
        sum(times(amplitudes, collectdnPl(order, cos(theta))), 1) + ...
        sin(theta)^3 .* sum(times(amplitudes, collectdnPl(order, cos(theta), 2)), 1);

    theta = pi/2 + pi/4;
    tol_theta = 1e-7;
    n = 1;
    % Newton Method!
    while abs(drdtheta(theta)) >= tol_theta && n < 150
        theta = mod(theta - drdtheta(theta)/dr2dtheta2(theta) - 1e-4, pi) + 1e-4; % If solution is close to pi, theta is unstable with mod function (therefore 1e-4 added)
        n = n + 1;
        if n == 50
            theta = 3.14159/2;
        elseif n == 100
            theta = rand()/100 + pi/2;
        end
        if n== 149
            println("Hey!")
        end
    end

    rmax = r_from_spherical(theta, amplitudes);

end