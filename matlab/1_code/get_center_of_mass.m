function [CM, contact_radius_guess, is_it_acceptable] = ...
        get_center_of_mass(amplitudes_guess, velocities_guess, previous_conditions, PROBLEM_CONSTANTS)

    n = length(previous_conditions); % Determines the order of the method
    if n > 2 || n < 1; throw("Hey!"); end

    extract_symbol = @(jj, field) previous_conditions{jj}.(field);
    amps_velocities = zeros(nb_harmonics, n+1);
    amps_deformations= zeros(nb_harmonics, n+1);
    for ii = 1:n
        amps_velocities(:, ii) = extract_symbol(ii, 'deformation_velocities');
        amps_deformations(:,ii) = extract_symbol(ii, 'deformation_amplitudes');
    end
    amps_deformations(:, end) = amplitudes_guess;
    amps_velocities(:, end) = velocities_guess;
    
    if n == 1
        coefs = [-1.0, 1.0];
    elseif n == 2
        rk = dt/previous_conditions{end}.dt;
        ak = (1+2*rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        coefs = [ck, bk, ak]; 
    end
    % There are three ways of calculating the new center of mass. 
    
    % Way 1: Using the relation z(t) = \sum_{l\geq 0} A_l P_l(pi)
    

    
end