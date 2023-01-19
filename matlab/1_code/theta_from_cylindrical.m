function tt = theta_from_cylindrical(r, amplitudes)
    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end

    zeta = zeta_generator(amplitudes);
    
    % Derivative of the function
    f_prime = @(theta) cos(theta) .* (1 + zeta(theta)) - sin(theta).^2 .* sum(times(amplitudes, collectdnPl(length(amplitudes), cos(theta))), 1);
    
    % Function to be minimized00
    f_objective = @(theta) sin(theta) * (1 + zeta(theta)) - r;

    theta = pi - 0.1;
    tol_theta = 1e-7;
    n = 1;

    % Newton Method!
    while abs(f_objective(theta)) >= tol_theta && n < 150
        theta = mod(theta - f_objective(theta)/f_prime(theta) - 1e-4, pi/2) + 1e-4 + pi/2; % If solution is close to pi, theta is unstable with mod function (therefore 1e-4 added)
        n = n + 1;
        if n == 50
            theta = 3.14159;
        elseif n == 100
            theta = rand() * pi/2 + pi/2;
        end
    end

    tt = theta;

end