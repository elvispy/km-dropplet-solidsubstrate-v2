function [probable_next_conditions, errortan] = ...
    get_next_step_v2(previous_conditions, dt, PROBLEM_CONSTANTS)
    % Hey!

    % Assume the tangency error is zero unless you find something wrong
    errortan = 0;
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
            errortan = error;
            return
        end
    end
    %is_it_acceptable = false;
    best_condition = {};
    minerror = inf;
    for idx = 1:100 
        [probable_next_conditions, is_it_acceptable, error] = ...
            advance_conditions(probable_next_conditions, ...
                previous_conditions, PROBLEM_CONSTANTS);
            
        if PROBLEM_CONSTANTS.DEBUG_FLAG
             plot_condition(2, probable_next_conditions, 1, ...
               struct("contact_radius", previous_conditions{end}.contact_radius, "iteration", idx));
        end
        if error < minerror
            minerror = error;
            best_condition = probable_next_conditions; 
        end

        if is_it_acceptable; break; end
        probable_next_conditions = ...
            pressure_with_flattened_surface(probable_next_conditions, ...
                previous_conditions, PROBLEM_CONSTANTS);
    end
    if idx == 100; probable_next_conditions = best_condition; end
    
    r_max = maximum_contact_radius(probable_next_conditions);
    %zeta = zeta_generator(probable_next_conditions);
    %z = @(theta) probable_next_conditions.center_of_mass +  cos(theta) .* (1 + zeta(theta));
    if r_max < probable_next_conditions.contact_radius || isnan(r_max)
        errortan = Inf;
        warning("Maximum contact radius passed!")
        return
    end
    
    if probable_next_conditions.contact_radius > 0
        rmax = dr * (probable_next_conditions.number_contact_points - 1/2);
        errortan = calculate_exit_angle(probable_next_conditions, ...
            theta_from_cylindrical(rmax, probable_next_conditions));
    elseif previous_conditions.contact_radius > PROBLEM_CONSTANTS.angle_tol 
        % If angle step is too big, we must make dt smaller
        errortan = inf;
    end
                
%         % Check that dropplet does not intersect with the substrate
%         for r_i = (dr * probable_next_conditions.number_contact_points):(dr/100):r_max
%             if z(theta_from_cylindrical(r_i, probable_next_conditions)) < -spatial_tol || isnan(r_i)
%                 errortan = Inf;
%                 break;
%             end
%         end
    if errortan < 0; disp("Errortan <= 0 ? "); end
end % end main function definition


function [new_probable_next_conditions, is_it_acceptable, error] = advance_conditions(probable_next_conditions, ...
    previous_conditions, PROBLEM_CONSTANTS)
    
    dt = probable_next_conditions.dt;
    
    [new_amplitudes, new_velocities] = solve_ODE_unkown(nan, probable_next_conditions, dt, previous_conditions, PROBLEM_CONSTANTS);
    
    nb_harmonics = previous_conditions{end}.nb_harmonics;
    
    n = length(previous_conditions); % Determines the order of the method
    if n > 2 || n < 1; throw("Hey!"); end

    extract_symbol = @(field) arrayfun(@(jj) previous_conditions{jj}.(field), 1:n);
    
    if n == 1
        coefs = [-1.0, 1.0];
    elseif n == 2
        rk = dt/previous_conditions{end}.dt;
        ak = (1+2*rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        coefs = [ck, bk, ak]; 
    end
    
    if norm(probable_next_conditions.pressure_amplitudes) == 0
        new_CM_velocity = (-dt/PROBLEM_CONSTANTS.froude_nb - sum(coefs(1:n) .* extract_symbol('center_of_mass_velocity')))/coefs(end);
        new_centerofmass = (dt * new_CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);    
    else
        % There are some ways to calculate the center of mass in case of contact. 
        switch PROBLEM_CONSTANTS.CM
            case 1 % First way: "Exact integration using pressure coefficients
                new_contact_radius = 0; % TODO: MODIFY THIS
                warning("This version is deprecated!");
                cos_theta = cos(theta_from_cylindrical(new_contact_radius, new_amplitudes)); %theta_from_cylindrical(previous_conditions{end}.contact_radius, new_amplitudes)); 
                d = 1 + sum(arrayfun(@(idx) (-1)^idx * new_amplitudes(idx), 2:nb_harmonics));
                d2zd2t = -1/PROBLEM_CONSTANTS.froude_nb - 3/2 *d^2 * ...
                    sum(probable_next_conditions.pressure_amplitudes .* ...
                        arrayfun(@(idx) integral(@(theta) my_legendre(idx, x)./x.^3, -1, cos_theta), 1:nb_harmonics));
                % (Now the B0 = -sum(Bl) term)
                f = @(x) -0.5 ./(x.^2);
                d2zd2t = d2zd2t -3/2 * d^2 * (-sum(probable_next_conditions.pressure_amplitudes) * (f(cos_theta) - f(-1)));
                new_CM_velocity  = (dt * d2zd2t - sum(coefs(1:n) .* extract_symbol('center_of_mass_velocity')))/coefs(end);
                new_centerofmass = (dt * new_CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);
            case 2 % Second way: Using the small oscillation's hypothesis
                warning("This version is deprecated because the first pressure amplitude is needed!");
                pressure_amplitudes_tentative = [probable_next_conditions.pressure_amplitudes, 0];
                Cl = @(l)  l * (l-1) / (2*l-1)     * pressure_amplitudes_tentative(l-1);
                Dl = @(l)  (l+2) * (l+1) / (2*l+3) * pressure_amplitudes_tentative(l+1);
                pressure_amplitudes_tentative = pressure_amplitudes_tentative(1:(end-1));

                d2zd2t = -1/PROBLEM_CONSTANTS.froude_nb - pressure_amplitudes_tentative(1);%- sum(coefs(1:n) .* extract_symbol("center_of_mass_velocity"), 1:n));

                for hb = 2:nb_harmonics
                    d2zd2t = d2zd2t + 3 * (new_amplitudes(hb) / (2*hb+1)) * (Cl(hb) - Dl(hb)); 
                end
                new_CM_velocity  = (dt * d2zd2t - sum(coefs(1:n) .* extract_symbol('center_of_mass_velocity')))/coefs(end);
                new_centerofmass = (dt * new_CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);
            case 3
                new_CM_velocity  = sum(arrayfun(@(idx) (-1)^idx * new_velocities(idx), 2:nb_harmonics));
                new_centerofmass = (dt * new_CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);
        end % end switch statement
    end
    
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
    
    my_tol = 1e-2; is_it_acceptable = false; 
    if norm(new_amplitudes - previous_conditions{end}.deformation_amplitudes) < my_tol ...
       && abs(new_contact_radius - previous_conditions{end}.contact_radius) < PROBLEM_CONSTANTS.angle_tol ...
       && error < my_tol
        is_it_acceptable = true;
    end
    
end


function [new_contact_radius, error] = calculate_contact_radius(new_amplitudes, ...
    new_centerofmass, PROBLEM_CONSTANTS, previous_contact_radius, probable_pressures)
    
    
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
    if sum(zs_evaluated < PROBLEM_CONSTANTS.spatial_tol)/N >= 0.995
        warning("Must make contact angle bigger!");
    end
    idx = 1;
    while zs_evaluated(idx) > PROBLEM_CONSTANTS.spatial_tol && idx < N
        idx = idx + 1;
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
    error = norm(zs_evaluated(min(idx, idx2):end));
    %else
    %    error = 
    %end
    %if any(intersect < spatial_tol); error = inf; end
        
end
% function new_probable_next_conditions = advance_conditions(probable_next_conditions, previous_conditions,  ...
%     new_nb_contact_points, dt, PROBLEM_CONSTANTS)
% 
%     % if previous_conditions <: ProblemConditions; previous_conditions = [previous_conditions]; end
%     n = length(previous_conditions); % Determines the order of the method
%     if n > 2 || n < 1; throw("Hey!"); end
%     extract_symbol = @(jj, field) previous_conditions{jj}.(field);
% 
%     nb_harmonics = previous_conditions{end}.nb_harmonics;
%     pressure_amplitudes_tentative = probable_next_conditions.pressure_amplitudes;
%     
%     if n == 1
%         coefs = [-1.0, 1.0];
%     elseif n == 2
%         rk = dt/previous_conditions{end}.dt;
%         ak = (1+2*rk)/(1+rk);
%         bk = -(1+rk);
%         ck = rk^2/(1+rk);
%         coefs = [ck, bk, ak]; 
%     end
% 
%     % Deformation amplitudes  at all necessary times
%     % dep_rhs = zeros(2, nb_harmonics);
%     amplitudes_tent = zeros(1, nb_harmonics);
%     amplitudes_velocities_tent = zeros(1, nb_harmonics);
%     
%     for hb = 2:nb_harmonics
%         dep_rhs = zeros(2, 1);
%         for jj = 1:n
%             dep_rhs = dep_rhs + coefs(jj) * ...
%                 [previous_conditions{jj}.deformation_amplitudes(hb); ...
%                  previous_conditions{jj}.deformation_velocities(hb)];
%         end
%         om = PROBLEM_CONSTANTS.omegas_frequencies(jj);
%         mat = [coefs(n+1), -dt; +dt*om^2, coefs(n+1)];
%         res = mat\(dt * [0; -hb * pressure_amplitudes_tentative(hb)] + ...
%                 dep_rhs);
%         amplitudes_tent(hb) = res(1);
%         amplitudes_velocities_tent(hb) = res(2);
%     end
% 
%     new_CM_velocity_times_dt = - sum(coefs(1:n) .* arrayfun(@(idx) extract_symbol(idx, "center_of_mass_velocity"), 1:n)) * dt;
%     
%     pressure_amplitudes_tentative = [pressure_amplitudes_tentative, 0];
%     Cl = @(l)  l * (l-1) / (2*l-1)     * pressure_amplitudes_tentative(l-1);
%     Dl = @(l)  (l+2) * (l+1) / (2*l+3) * pressure_amplitudes_tentative(l+1);
%     pressure_amplitudes_tentative = pressure_amplitudes_tentative(1:(end-1));
%     
%     new_CM_velocity_times_dt = new_CM_velocity_times_dt - dt^2 * PROBLEM_CONSTANTS.froude_nb; 
%     %Special case: First harmonics
%     new_CM_velocity_times_dt = new_CM_velocity_times_dt - dt^2 * pressure_amplitudes_tentative(1); 
%     % General case
%     
%     for hb = 2:nb_harmonics
%         new_CM_velocity_times_dt = new_CM_velocity_times_dt + 3 * dt^2 * ...
%             (amplitudes_tent(hb) / (2*hb+1)) * (Cl(hb) - Dl(hb)); 
%     end
%     %new_CM_velocity_times_dt = new_CM_velocity_times_dt + 3 * dt^2 * ...
%     %        (amplitudes_tent(ll) / (2*ll+1)) * (Cl(ll) - 0); 
%     
% 
%     new_CM_velocity_times_dt = new_CM_velocity_times_dt / coefs(end);
%     
%     new_center_of_mass = (new_CM_velocity_times_dt - sum(coefs(1:n) .* arrayfun(@(idx) extract_symbol(idx, "center_of_mass"), 1:n)));
%     new_center_of_mass = new_center_of_mass /coefs(end);
% 
%     new_probable_next_conditions = ProblemConditions( ...
%         nb_harmonics,...
%         amplitudes_tent,...
%         amplitudes_velocities_tent,...
%         pressure_amplitudes_tentative,...
%         previous_conditions{end}.current_time + dt,...
%         dt,...
%         new_center_of_mass,...
%         new_CM_velocity_times_dt/dt, ...% Divide by dt!
%         new_nb_contact_points);
% 
% end



function [probable_next_conditions, is_it_acceptable, previous_tentatives, idxs] = ...
    update_tentative_heuristic(probable_next_conditions, ~, ... % previous_conditions in place of ~
    dr, ~, spatial_tol, PROBLEM_CONSTANTS, previous_tentatives, idxs) % dt in place of ~
    
    harmonics_qtt = probable_next_conditions.nb_harmonics;
    NB_SAMPLES = harmonics_qtt + 1;
    is_it_acceptable = true;
    % First, lets check if the given probable next condition is acceptable or not.
    heights = probable_next_conditions.center_of_mass * ones(1, NB_SAMPLES);
    pressure_samples = zeros(1, NB_SAMPLES);
    pressure_amps = zeta_generator(probable_next_conditions.pressure_amplitudes);

    zeta = zeta_generator(probable_next_conditions);
    contact_radius = dr * (probable_next_conditions.number_contact_points - 1/2);

    % Check if heights are in bounds
    for ii = 1:NB_SAMPLES
        % Angle of last contact point
        theta = theta_from_cylindrical(contact_radius*(ii-1)/(NB_SAMPLES-1), probable_next_conditions.deformation_amplitudes);
        heights(ii) = heights(ii) +  cos(theta) * (1 + zeta(theta)); 
        pressure_samples(ii) = pressure_amps(theta) - sum(probable_next_conditions.pressure_amplitudes);
        if abs(heights(ii)) > spatial_tol
            is_it_acceptable = false; % We dont break because we will need all heights to guess a new pressure profile
            % if PROBLEM_CONSTANTS.DEBUG_FLAG; disp("Breaking because heights do not conform to tolerances"); end
        end
    end
    
    % We dont care about last point
    %pressure_samples(end) = 0;
    %heights(end) = 0;
    
    % you cant have negative pressures
    % assert( min(pressure_samples) >= 0,  "Negative pressures!");
    err = norm(heights);
    added = false;
    if length(previous_tentatives) < 2 || err < previous_tentatives{end}.error || err < previous_tentatives{end-1}.error
        previous_tentatives = {previous_tentatives{1:end} struct("heights", heights,"pressure_samples", pressure_samples, "error", err)};
        
        n = length(previous_tentatives);
        if n > 2
            added = true;
            if err < previous_tentatives{idxs(1)}.error && previous_tentatives{idxs(1)}.error > previous_tentatives{idxs(2)}.error
                idxs(1) = length(previous_tentatives);
            elseif err < previous_tentatives{idxs(2)}.error
                idxs(2) = length(previous_tentatives);    
            end
        end
    end
        
    % If some height is out of bound, try to adapt pressure coefficients
    if is_it_acceptable == false
        % Heuristic tentative: Increase of reduce the pressure at given points to fit flat area.
        theta_max = theta_from_cylindrical(contact_radius, probable_next_conditions.deformation_amplitudes);
            
        % y_velocity(theta::Float64)   = cos(theta)^2 * sum(probable_next_conditions.deformation_velocities .* 
        %         (collectdnPl(cos(theta); lmax = order, n = 1).parent));
        %r_positions = ?r * (0:(probable_next_conditions.new_number_contact_points-1));

        % Perturbation at LinRange(0, rmax, harmonics_qtt) to flatten the surface
        % pressure_perturbation = zeros(harmonics_qtt, 1);

        if length(previous_tentatives) < 2
            % Modify pressure amplitudes so as to try to flatten the surface
            pressure_perturbation = times(((heights >= 0) - 0.5), -0.2 * abs(pressure_samples)) ...%arrayfun(@(h) 0.1 * (h+eps)/abs(h+eps), heights)) + ...
                + (heights <-spatial_tol) .* (abs(pressure_samples) < 1e-3/PROBLEM_CONSTANTS.pressure_unit) *  .05/PROBLEM_CONSTANTS.pressure_unit ...
                + (heights > spatial_tol) .* (abs(pressure_samples) < 1e-3/PROBLEM_CONSTANTS.pressure_unit) * -.01/PROBLEM_CONSTANTS.pressure_unit ; %arrayfun(@(idx) (heights(idx) < 0) *  ( abs(pressure_samples(idx)) < 1e-8) * 0.01, 1:harmonics_qtt);
            %theta = theta_from_cylindrical(rmax*(ii-1)/(harmonics_qtt-1), probable_next_conditions.deformation_amplitudes)
            %%-pressure_perturbation(end) = 0;
        else
            % First tentative: assume linearity between the last two tentatives
            % This interpolator tries to 
            
            if added == false %%&& rand() > 0.4
                pressure_perturbation = times(-abs(previous_tentatives{idxs(2)}.pressure_samples), rand()/15 * ((heights >= 0) - 0.5));
            else
                interpolator = @(d1, d2, idx)  ... 
                     ((d2.heights(idx) * d1.pressure_samples(idx) ...
                    - d1.heights(idx) * d2.pressure_samples(idx)) / ...
                     (d2.heights(idx) - d1.heights(idx)));
                d1 = previous_tentatives{idxs(1)};
                d2 = previous_tentatives{idxs(2)};

                % We only apply perturbe pressure where the heights is
                % unnacceptable, 
                pressure_perturbation = (abs(heights) > spatial_tol/5) .* ...
                    (arrayfun(@(idx) interpolator(d1, d2, idx), 1:NB_SAMPLES) - pressure_samples);
                idxs_2 = and(heights < -spatial_tol, pressure_perturbation < 0);
                pressure_perturbation(idxs_2) = abs(pressure_samples(idxs_2)) * rand()/15;
                %%-pressure_perturbation(end) = 0;
            end
        end

        % Compute new pressure coefficients:

        % Interpolate linearly between pressure points

        f = @(r) interp1(contact_radius * linspace(0, 1, NB_SAMPLES), pressure_perturbation, r, 'linear',  0); 
        ps = @(theta) f(r_from_spherical(theta, probable_next_conditions.deformation_amplitudes));

        projected_pressure_perturbations = project_amplitudes(ps, harmonics_qtt, [theta_max, pi], PROBLEM_CONSTANTS, true);
        
        % Now let's kill pressure distribution outside of contact area
        %%-pressure_outside_perturbation = kill_pressure_outside(probable_next_conditions, theta_max);

        %assert(all(probable_next_conditions.pressure_amplitudes + projected_pressure_amplitudes >= 0),  "Need to fix this");

        probable_next_conditions.pressure_amplitudes = probable_next_conditions.pressure_amplitudes ...
            + projected_pressure_perturbations;
%         probable_next_conditions = ProblemConditions( ...
%             probable_next_conditions.nb_harmonics, ...
%             probable_next_conditions.deformation_amplitudes, ...
%             probable_next_conditions.deformation_velocities, ...
%             probable_next_conditions.pressure_amplitudes + projected_pressure_perturbations, ...
%             probable_next_conditions.current_time, ...
%             probable_next_conditions.dt, ...
%             probable_next_conditions.center_of_mass, ...
%             probable_next_conditions.center_of_mass_velocity, ...
%             probable_next_conditions.number_contact_points);
   
    end
    
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


