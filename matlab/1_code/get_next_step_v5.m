function [probable_next_conditions, errorflag] = ...
    get_next_step_v5(previous_conditions, dt, PROBLEM_CONSTANTS)
    % This method tries to minimize the objective function via Differential
    % Evoution
    maxtry = 5;
    Flat = ceil(previous_conditions{end}.nb_harmonics*0.9);
    Ps = previous_conditions{end}.nb_harmonics - Flat;
    ftm = @(Xn, theta_contact) function_to_minimize(Xn, previous_conditions, dt, theta_contact, PROBLEM_CONSTANTS.froude_nb, Ps, Flat);
    jaccalc = @(Xn, theta_contact) JacobianCalculator(Xn, previous_conditions, dt, theta_contact, Ps, Flat);
    % Build initial conditions
    A = NaN * ones(maxtry, 2);
    conditions = cell(1, maxtry);
    
    
    for current_guess_idx = -2:2 %1:maxtry
        theta_contact = previous_conditions{end}.contact_angle + PROBLEM_CONSTANTS.angle_tol * current_guess_idx;
        if theta_contact > pi
            if A(current_guess_idx+2, 1) == pi
                A(current_guess_idx+3, :) = A(current_guess_idx+2, :);
                conditions{current_guess_idx+3} = conditions{current_guess_idx+2};
                continue
            else            
                theta_contact = pi;
            end
        end
        
        % Newton Method
        Xn = [previous_conditions{end}.deformation_amplitudes(2:end)'; previous_conditions{end}.deformation_velocities(2:end)'; ...
                previous_conditions{end}.pressure_amplitudes(:); previous_conditions{end}.center_of_mass; ...
                previous_conditions{end}.center_of_mass_velocity]; % Initial guess for Newton-Raphson method
        Xmin = Xn;
        X1 = Xn;
        fmin = Inf;
        EPS = 10;
        for m = 1:100

            valvec = ftm(Xn, theta_contact);
            normvalvec = norm(valvec);
               
            if normvalvec/fmin > 1e+5 
                % Restart or redo
                break
            elseif normvalvec < 1e-3
                Xmin = Xn;
                fmin = normvalvec;
                break
            elseif normvalvec < fmin
                Xmin = Xn;
                fmin = normvalvec;
            end
            
            Jac = jaccalc(Xn, theta_contact);
            %myres = lsqr(Jac, valvec, 1e-3, 100, [], [], valvec);
            %Xn = Xn - myres;
            % If solution has gone too far, try to stabilize by choosing
            % one near.
            nn = norm(Xn-X1);
            if nn > EPS
               Xn = X1 + (Xn-X1)/nn * rand();
            end
        end
        
        A(current_guess_idx+3, :) = [theta_contact, fmin];
        conditions{current_guess_idx+3} = Xmin;
        
        
        % Update theta_contact

    end
    [~, idx] = min(A(:, 2));
    Xn = conditions{idx};
    if abs(idx-3) > 1
        errorflag = true;
    else
        errorflag = false;
    end
    % Lets try with the same pressure
    probable_next_conditions = previous_conditions{end};
    N = probable_next_conditions.nb_harmonics;
    probable_next_conditions.deformation_amplitudes = [0; Xn(1:(N-1))]';
    probable_next_conditions.deformation_velocities = [0; Xn(N:(2*N-2))]';
    probable_next_conditions.pressure_amplitudes    = Xn((2*N-1):(end-2))';
    probable_next_conditions.center_of_mass         = Xn(end-1);
    probable_next_conditions.center_of_mass_velocity= Xn(end);
    probable_next_conditions.dt = dt;
    probable_next_conditions.current_time = previous_conditions{end}.current_time + dt;
    probable_next_conditions.contact_angle = A(idx, 1);
    probable_next_conditions.contact_radius = round(sin(A(idx, 1)) * (1 + ...
        sum(arrayfun(@(idx) (-1)^(idx+1) * Xn(idx), 1:length(Xn)))), 10);

    
end % end main function definition

function contact_angle = calculate_contact_angle(new_amplitudes, ...
    new_centerofmass, PROBLEM_CONSTANTS)%, previous_contact_radius, probable_pressures)
    

    
    zeta = zeta_generator(new_amplitudes);
    z = @(theta) new_centerofmass + cos(theta) .* (1 + zeta(theta));
    
    N = 500;
    contact_angle_upper = pi;
    contact_angle_lower = pi/2;
    while contact_angle_upper - contact_angle_lower >= PROBLEM_CONSTANTS.angle_tol
        midpoint = contact_angle_upper/2 + contact_angle_lower/2; 
        if z(midpoint) < 0
            contact_angle_upper = midpoint;
        else
            contact_angle_lower = midpoint;
        end
    end
    if all(z(( pi - contact_angle_upper ) * rand(20, 1) + contact_angle_upper) < 0)
        contact_angle = contact_angle_upper;
    else
        throw("maximum contact angle not found")
    end
    
        
end