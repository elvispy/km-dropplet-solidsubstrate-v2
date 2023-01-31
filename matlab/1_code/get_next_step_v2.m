function [probable_next_conditions, errorflag] = ...
    get_next_step_v2(previous_conditions, dt, PROBLEM_CONSTANTS)

    % Assume the tangency error is zero unless you find something wrong
    errorflag = false;
    %spatial_tol = PROBLEM_CONSTANTS.spatial_tol;
    %angle_tol   = PROBLEM_CONSTANTS.angle_tol;
    % Lets try with the same pressure
    probable_next_conditions = previous_conditions{end};
    probable_next_conditions.current_time = probable_next_conditions.current_time + dt;
    probable_next_conditions.dt = dt;
    if probable_next_conditions.center_of_mass_velocity > 0 ...% If there is chance of not intersecting in next time
         && previous_conditions{end}.contact_radius > pi - PROBLEM_CONSTANTS.angle_tol
        no_pressure_conditions = probable_next_conditions;
        no_pressure_conditions.pressure_amplitudes = zeros(probable_next_conditions.nb_harmonics, 1);
        [no_pressure_conditions, ~, error] = advance_conditions(no_pressure_conditions, ...
            previous_conditions, PROBLEM_CONSTANTS);
        
        if error == 0
            probable_next_conditions = no_pressure_conditions;
            errorflag = false;
            return
        end
    end
    %is_it_acceptable = false;
    best_condition = {};
    minerror = inf;
    for idx = 1:50
        [probable_next_conditions, is_it_acceptable, error] = ...
            advance_conditions(probable_next_conditions, ...
                previous_conditions, PROBLEM_CONSTANTS);
            
        if PROBLEM_CONSTANTS.DEBUG_FLAG
             plot_condition(2, probable_next_conditions, 1, ...
               struct("contact_radius", previous_conditions{end}.contact_radius, "iteration", idx));
             %sp = zeta_generator(probable_next_conditions);
             %fprintf("dt:%.5g, south pole: %.5g\n", probable_next_conditions.dt, ...
             %    probable_next_conditions.center_of_mass - (1 + sp(pi)));
        end
        if error < minerror
            minerror = error;
            best_condition = probable_next_conditions; 
        end

        if is_it_acceptable; break; end
        probable_next_conditions = ...
            update_pressure_guess(probable_next_conditions, ...
                previous_conditions, PROBLEM_CONSTANTS);
    end
    if idx == 50; probable_next_conditions = best_condition; errorflag = true; end
    
    r_max = maximum_contact_radius(probable_next_conditions);
    %zeta = zeta_generator(probable_next_conditions);
    %z = @(theta) probable_next_conditions.center_of_mass +  cos(theta) .* (1 + zeta(theta));
    if r_max < probable_next_conditions.contact_radius || isnan(r_max)
        errorflag = true;
        warning("Maximum contact radius passed!")
        return
    end
    
    if probable_next_conditions.contact_radius < 1e-6 && norm(probable_next_conditions.pressure_amplitudes) > 10
        errorflag = true;
    end
    
%     if probable_next_conditions.contact_radius > 0
%         % rmax = maximum_contact_radius(probable_next_conditions); % dr * (probable_next_conditions.number_contact_points - 1/2); %% TODO FIX THIS
%         % errorflag = (abs(best_condition.pressure_deformation(1)) > 1e+4);
%     else
    if abs(theta_from_cylindrical( probable_next_conditions.contact_radius, probable_next_conditions)...
           - theta_from_cylindrical(previous_conditions{end}.contact_radius, previous_conditions{end})) > PROBLEM_CONSTANTS.angle_tol
        % If angle step is too big, we must make dt smaller
        errorflag = true;
    end
                
%         % Check that dropplet does not intersect with the substrate
%         for r_i = (dr * probable_next_conditions.number_contact_points):(dr/100):r_max
%             if z(theta_from_cylindrical(r_i, probable_next_conditions)) < -spatial_tol || isnan(r_i)
%                 errortan = Inf;
%                 break;
%             end
%         end
    
end % end main function definition


function [new_probable_next_conditions, is_it_acceptable, error] = advance_conditions(probable_next_conditions, ...
    previous_conditions, PROBLEM_CONSTANTS)
    
    dt = probable_next_conditions.dt;
    
    [new_amplitudes, new_velocities] = solve_ODE_unkown(nan, probable_next_conditions, dt, previous_conditions, PROBLEM_CONSTANTS);
    
    nb_harmonics = previous_conditions{end}.nb_harmonics;
    
    [new_centerofmass, new_CM_velocity] = get_center_of_mass(new_amplitudes, new_velocities, ...
        probable_next_conditions, previous_conditions, PROBLEM_CONSTANTS);

    [new_contact_radius, error] = calculate_contact_radius(new_amplitudes, new_centerofmass, ...
        PROBLEM_CONSTANTS, previous_conditions{end}.contact_radius, probable_next_conditions.pressure_amplitudes);
        
    new_probable_next_conditions = ProblemConditions_v2( ...
    nb_harmonics,...
    new_amplitudes,...
    new_velocities,...
    probable_next_conditions.pressure_amplitudes,...
    previous_conditions{end}.current_time + dt,...
    dt,...
    new_centerofmass,...
    new_CM_velocity, ...
    new_contact_radius);

    % Now let's check if convergence was attained
    my_tol = 1e-4; is_it_acceptable = false; 
    if norm(new_amplitudes - probable_next_conditions.deformation_amplitudes) < my_tol ... %&& abs(new_contact_radius - probable_next_conditions.contact_radius) < my_tol ...%r_from_spherical(PROBLEM_CONSTANTS.angle_tol, probable_next_conditions) ...
       && error < my_tol
        sp = zeta_generator(probable_next_conditions);
        if ~(probable_next_conditions.center_of_mass - (1 + sp(pi)) > PROBLEM_CONSTANTS.spatial_tol ...
                && new_contact_radius < 1e-6)
            is_it_acceptable = true;
        end
    end
end % end advance conditioins definition


function [new_contact_radius, error] = calculate_contact_radius(new_amplitudes, ...
    new_centerofmass, PROBLEM_CONSTANTS, previous_contact_radius, probable_pressures)
    
%     angle = pi;
%     crossed = false;
%     step = pi/500;
%     tol = pi/100000;
%     down = true;
%     while step >= tol
%         val = z(angle);
%         if val < -PROBLEM_CONSTANTS.spatial_tol; crossed = true; end
%         
%         if down == true && val >= -PROBLEM_CONSTANTS.spatial_tol
%             down = false;
%             step = step / 2;
%         elseif down == false && val < -PROBLEM_CONSTANTS.spatial_tol
%             down = true;
%             step = step / 2;
%         end
%         if down == true
%             angle = angle - step; 
%         else
%             angle = angle + step;
%         end
%     end    
%     
%     if crossed == true && norm(probable_pressures) == 0
%         error = inf;
%     elseif norm(probable_pressures) == 0
%         error = 0;
%     end
    
    zeta = zeta_generator(new_amplitudes);
    z = @(theta) new_centerofmass + cos(theta) .* (1 + zeta(theta));
    
    N = 500;
    thetas = linspace(max(pi/2, theta_from_cylindrical(previous_contact_radius, new_amplitudes) - 2*pi/10), pi, N);
    
    zs_evaluated = z(thetas);    
    
    if norm(probable_pressures) == 0
        
        error = 1e+10 * any(zs_evaluated < - PROBLEM_CONSTANTS.spatial_tol);
        if error == 0
            new_contact_radius = 0;
            return
        end
    end
    % Penalizing places with too little pressure
    if sum(zs_evaluated < 0.0*PROBLEM_CONSTANTS.spatial_tol)/N >= 0.995
        warning("Must make contact angle bigger!");
    end
    idx = N;
    while zs_evaluated(idx) <= 0.0 * PROBLEM_CONSTANTS.spatial_tol && idx > 1
        idx = idx - 1;
    end
    
    % Penalizing places where there is pressure but no contact
%     Ps = zeta_generator(probable_pressures);
%     ps_evaluated = abs(Ps(thetas) - sum(probable_pressures));
%     ps_evaluated = (ps_evaluated+eps)/(max(ps_evaluated)+eps); % Normalizing pressure
%     
%     if  sum(ps_evaluated > 0.01)/N > 0.995 && norm(probable_pressures) > 0
%         warning("Must make contact angle bigger!");
%     end
%     if min(ps_evaluated) > 1/20 % If the ratio min/max (max = 1) is too big, means that there may be no contact.
%         idx2 = N;
%     else
%         idx2 = 1;
%         while ps_evaluated(idx) < 0.2 && idx2 < N
%             idx2 = idx2 + 1;
%         end
%     end
    idx2 = N;
    
    %if sum(abs(intersect(idx:end)) < spatial_tol)/N >= 0.995
    new_contact_radius = r_from_spherical(thetas(min(idx, idx2)), new_amplitudes);
    %%-crp = r_from_spherical(theta_from_cylindrical(previous_contact_radius, new_amplitudes)-PROBLEM_CONSTANTS.angle_tol, new_amplitudes);
    %%-new_contact_radius = min(crp, new_contact_radius);
    error = norm(zs_evaluated(min(idx, idx2):end));

    if pi - thetas(idx) <= PROBLEM_CONSTANTS.angle_tol
        new_contact_radius = 0.0;
    end
    %else
    %    error = 
    %end
    %if any(intersect < spatial_tol); error = inf; end
        
end




%Calculates the exit angle at the contact angle theta

function exit_angle = calculate_exit_angle(amplitudes, angle)
    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    
    zeta = zeta_generator(amplitudes);
    if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
    der = @(theta) sum(amplitudes .* collectdnPl(length(amplitudes), cos(theta)), 1);

    dzdr = @(theta)  (-sin(theta) .* (1 + zeta(theta)) - cos(theta) .* sin(theta) .* der(theta)) ./ ...
        (cos(theta) .* (1 + zeta(theta)) - sin(theta).^2 .* der(theta));
        
    exit_angle =  dzdr(angle);
end


