function new_probable_next_conditions = update_pressure_guess(probable_next_conditions, ...
            previous_conditions, PROBLEM_CONSTANTS)
        
    if probable_next_conditions.contact_radius < 1e-6
        pressure_guess = (0.5 + rand()/2) * probable_next_conditions.pressure_amplitudes;
        new_probable_next_conditions = probable_next_conditions;
        new_probable_next_conditions.pressure_amplitudes = pressure_guess;
        return
    end
        
    theta_contact = theta_from_cylindrical(probable_next_conditions.contact_radius, probable_next_conditions);
    
    switch PROBLEM_CONSTANTS.PG
        case 1
            % Least squares on flat surface B(1) is ignored
            
            nb_harmonics = probable_next_conditions.nb_harmonics;
            M = 2 * nb_harmonics;
            thetas = linspace(theta_contact, pi, M);
            LVAL = collectPl(nb_harmonics, cos(thetas));
            LVAL = LVAL(2:end, :)';
            b = -probable_next_conditions.center_of_mass./cos(thetas) - 1; b = b';
            [amps, ~] = lsqr(LVAL, b, [], 100, [], [], probable_next_conditions.deformation_amplitudes(2:end)');
            amplitudes_modified = [0, amps'];

        case 2
            % Flatten by "Frankenstein"
            nb_harmonics = probable_next_conditions.nb_harmonics;

            %theta_contact = theta_from_cylindrical(probable_next_conditions.contact_radius, probable_next_conditions);

            g = zeta_generator(probable_next_conditions);
            f = @(theta) f_generator(theta, theta_contact, g, probable_next_conditions.center_of_mass);
            amplitudes_modified = project_amplitudes(f, nb_harmonics, [0, pi, theta_contact], PROBLEM_CONSTANTS, true);
            amplitudes_modified(1) = 0;   

        case 3
            % Least squares on flat surface B(1) is ignored
            
            nb_harmonics = probable_next_conditions.nb_harmonics;
            M = 2 * nb_harmonics;
            thetas = logspace(log10(theta_contact), log10(pi), M); % linspace(theta_contact, pi, M);
            LVAL = collectPl(nb_harmonics, cos(thetas));
            LVAL = LVAL(2:end, :)';
            b = -probable_next_conditions.center_of_mass./cos(thetas) - 1; b = b';
            [amps, ~] = lsqr(LVAL, b, [], 100, [], [], probable_next_conditions.deformation_amplitudes(2:end)');
            amplitudes_modified = [0, amps'];

    end
    
    pressure_guess = solve_ODE_unkown(amplitudes_modified, nan, probable_next_conditions.dt, previous_conditions, PROBLEM_CONSTANTS);
    
    new_probable_next_conditions = probable_next_conditions;
    new_probable_next_conditions.pressure_amplitudes = pressure_guess;
    if PROBLEM_CONSTANTS.KILL_OUTSIDE == true
        new_probable_next_conditions.pressure_amplitudes = kill_pressure_outside_slow(new_probable_next_conditions, theta_contact, PROBLEM_CONSTANTS);
    end
end

function v = f_generator(theta, theta_contact, zeta, d)
    v = (theta <= theta_contact) .* (1 + zeta(theta)) ...
      + (theta >  theta_contact) .* -d./cos(theta); %+ zeros(size(theta));
    
%     for ii = 1:length(theta)
%         if theta(ii) <= theta_contact
%             v(ii) = 1 + zeta(theta(ii));
%         else
%             v(ii) = -d./cos(theta(ii));
%         end
%     end
end