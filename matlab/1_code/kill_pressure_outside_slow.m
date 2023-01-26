function pressure_amplitudes = kill_pressure_outside_slow(probable_next_conditions, theta_max, ~)
    
    N = probable_next_conditions.nb_harmonics;

    %Past pressure coeficients 
    Bl = probable_next_conditions.pressure_amplitudes;
    
    % Pressure perturbation coefficients
    %Cl = zeros(1, N);
    z = zeta_generator(Bl);
    %avg = mean(z(linspace(0, pi/10)));
    %f = @(ang) (ang > theta_max) .* (z(ang) - avg);
    
    pressure_amplitudes = project_amplitudes(z, N, [theta_max, pi], true, true);
    if size(pressure_amplitudes, 1) > 1; pressure_amplitudes = pressure_amplitudes'; end
end