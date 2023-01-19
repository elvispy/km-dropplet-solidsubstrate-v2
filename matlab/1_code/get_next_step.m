function [probable_next_conditions, errortan] = ...
    get_next_step(previous_conditions, new_number_contact_points, dt, ...
        dr, spatial_tol, PROBLEM_CONSTANTS)
    

    if new_number_contact_points < 0
        errortan = Inf;
        probable_next_conditions = struct( ...
        "harmonics_qtt", NaN, "deformation_amplitudes", NaN, ...
        "deformation_velocities", NaN, "pressure_deformation", NaN, ...
        "current_time", NaN, "dt", NaN, "center_of_mass", NaN, "center_of_mass_velocity", NaN, "number_contact_points", NaN);
        return
        
    else
        % Assume the tangency error is zero unless you find something wrong
        errortan = 0;
        % Lets try with no distribution
        probable_next_conditions = ProblemConditions( ...
                NaN, NaN, NaN, zeros(1, previous_conditions{end}.nb_harmonics), NaN, NaN, NaN, NaN, NaN);
        probable_next_conditions = advance_conditions(probable_next_conditions, ...
                previous_conditions, new_number_contact_points, dt, PROBLEM_CONSTANTS);
            
        if PROBLEM_CONSTANTS.DEBUG_FLAG
           plot_condition(2, probable_next_conditions, dr, 3, ...
               struct("previous_contact_points", previous_conditions{end}.number_contact_points, "iteration", 1));
        end
            
        f = zeta_generator(probable_next_conditions);
        thetas = arrayfun(@(r) theta_from_cylindrical(r, probable_next_conditions), ...
            linspace(0, (new_number_contact_points  - 1/2) * dr, probable_next_conditions.nb_harmonics)); 
        % If with no pressure we still dont descend enough, this number of
        % contact points is no good
        heights = probable_next_conditions.center_of_mass + cos(thetas) .* (1 + f(thetas));
                % Try with a similar distribution
                %         if new_number_contact_points > 0
                %             probable_next_conditions = ProblemConditions( ...
                %                 NaN, NaN, NaN, previous_conditions{end}.pressure_amplitudes, NaN, NaN, NaN, NaN, NaN);
                %         else
                %             probable_next_conditions = ProblemConditions( ...
                %                 NaN, NaN, NaN, zeros(1, previous_conditions{end}.nb_harmonics), NaN, NaN, NaN, NaN, NaN);
                %         end
        
        if any(heights > 1/50) && new_number_contact_points > 0 %%-&& false %% TODO REVERT THIS
            [probable_next_conditions, errortan] = get_next_step(probable_next_conditions, -1, NaN, NaN, NaN, NaN);
            return
        elseif new_number_contact_points > 0
            % Save previous tentative
            if norm(probable_next_conditions.pressure_amplitudes) > 0
                [~, ~, previous_tentatives, ~] = ...
                        update_tentative_heuristic(probable_next_conditions, previous_conditions, ...
                            dr, dt, spatial_tol, PROBLEM_CONSTANTS, {}, []);
            else
                previous_tentatives = {};
            end
            % Try with same pressure distribution
            
            probable_next_conditions = ProblemConditions( ...
                NaN, NaN, NaN, previous_conditions{end}.pressure_amplitudes, NaN, NaN, NaN, NaN, NaN);
            
            iteration = 0;
            idxs = [1, 2];
            %previous_tentatives = {};
            %is_it_acceptable = true;
            while iteration < 50  % 100
                iteration = iteration + 1;

                probable_next_conditions = advance_conditions(probable_next_conditions, ...
                    previous_conditions, new_number_contact_points, dt, PROBLEM_CONSTANTS);
                
                if PROBLEM_CONSTANTS.DEBUG_FLAG
                   plot_condition(2, probable_next_conditions, dr, 3, ...
                       struct("previous_contact_points", previous_conditions{end}.number_contact_points, "iteration", iteration));
                end
                %if new_number_contact_points > 0 % Only update pressures if you expect to have pressure
                        [probable_next_conditions, is_it_acceptable, previous_tentatives, idxs] = ...
                        update_tentative_heuristic(probable_next_conditions, previous_conditions, ...
                            dr, dt, spatial_tol, PROBLEM_CONSTANTS, previous_tentatives, idxs);
                %end

                if is_it_acceptable == true; break;  end
                if iteration == 50; disp("Hey!, we reached 50 iterations."); end
            end

        end

        %if new_number_contact_points == 0; errortan = 0; end
        % Check if solution makes sense
        r_max = maximum_contact_radius(probable_next_conditions);
        zeta = zeta_generator(probable_next_conditions);
        z = @(theta) probable_next_conditions.center_of_mass +  cos(theta) .* (1 + zeta(theta));
        if r_max < dr * (probable_next_conditions.number_contact_points - 1/2) || isnan(r_max)
            [probable_next_conditions, errortan] = get_next_step(probable_next_conditions, -1, NaN, NaN, NaN, NaN);
        end

        % Calculate errortan
        if new_number_contact_points > 0
            rmax = dr * (probable_next_conditions.number_contact_points - 1/2);
            errortan = calculate_exit_angle(probable_next_conditions, ...
                theta_from_cylindrical(rmax, probable_next_conditions));
        end
        
        % Check that dropplet does not intersect with the substrate
        for r_i = (dr * probable_next_conditions.number_contact_points):(dr/100):r_max
            if z(theta_from_cylindrical(r_i, probable_next_conditions)) < -spatial_tol || isnan(r_i)
                errortan = Inf;
                break;
            end
        end

    end % End main if
    if errortan < 0; disp("Errortan <= 0 ? "); end
end % end main function definition


function new_probable_next_conditions = advance_conditions(probable_next_conditions, previous_conditions,  ...
    new_nb_contact_points, dt, PROBLEM_CONSTANTS)

    % if previous_conditions <: ProblemConditions; previous_conditions = [previous_conditions]; end
    n = length(previous_conditions); % Determines the order of the method
    if n > 2 || n < 1; throw("Hey!"); end
    extract_symbol = @(jj, field) previous_conditions{jj}.(field);

    nb_harmonics = previous_conditions{end}.nb_harmonics;
    pressure_amplitudes_tentative = probable_next_conditions.pressure_amplitudes;
    
    if n == 1
        coefs = [-1.0, 1.0];
    elseif n == 2
        rk = dt/previous_conditions{end}.dt;
        ak = (1+2*rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        coefs = [ck, bk, ak]; 
    end

    % Deformation amplitudes  at all necessary times
    % dep_rhs = zeros(2, nb_harmonics);
    amplitudes_tent = zeros(1, nb_harmonics);
    amplitudes_velocities_tent = zeros(1, nb_harmonics);
    
    for hb = 2:nb_harmonics
        dep_rhs = zeros(2, 1);
        for jj = 1:n
            dep_rhs = dep_rhs + coefs(jj) * ...
                [previous_conditions{jj}.deformation_amplitudes(hb); ...
                 previous_conditions{jj}.deformation_velocities(hb)];
        end
        om = PROBLEM_CONSTANTS.omegas_frequencies(jj);
        mat = [coefs(n+1), -dt; +dt*om^2, coefs(n+1)];
        res = mat\(dt * [0; -hb * pressure_amplitudes_tentative(hb)] + ...
                dep_rhs);
        amplitudes_tent(hb) = res(1);
        amplitudes_velocities_tent(hb) = res(2);
    end

    new_CM_velocity_times_dt = - sum(coefs(1:n) .* arrayfun(@(idx) extract_symbol(idx, "center_of_mass_velocity"), 1:n)) * dt;
    
    pressure_amplitudes_tentative = [pressure_amplitudes_tentative, 0];
    Cl = @(l)  l * (l-1) / (2*l-1)     * pressure_amplitudes_tentative(l-1);
    Dl = @(l)  (l+2) * (l+1) / (2*l+3) * pressure_amplitudes_tentative(l+1);
    pressure_amplitudes_tentative = pressure_amplitudes_tentative(1:(end-1));
    
    new_CM_velocity_times_dt = new_CM_velocity_times_dt - dt^2 * PROBLEM_CONSTANTS.froude_nb; 
    %Special case: First harmonics
    new_CM_velocity_times_dt = new_CM_velocity_times_dt - dt^2 * pressure_amplitudes_tentative(1); 
    % General case
    
    for hb = 2:nb_harmonics
        new_CM_velocity_times_dt = new_CM_velocity_times_dt + 3 * dt^2 * ...
            (amplitudes_tent(hb) / (2*hb+1)) * (Cl(hb) - Dl(hb)); 
    end
    %new_CM_velocity_times_dt = new_CM_velocity_times_dt + 3 * dt^2 * ...
    %        (amplitudes_tent(ll) / (2*ll+1)) * (Cl(ll) - 0); 
    

    new_CM_velocity_times_dt = new_CM_velocity_times_dt / coefs(end);
    
    new_center_of_mass = (new_CM_velocity_times_dt - sum(coefs(1:n) .* arrayfun(@(idx) extract_symbol(idx, "center_of_mass"), 1:n)));
    new_center_of_mass = new_center_of_mass /coefs(end);

    new_probable_next_conditions = ProblemConditions( ...
        nb_harmonics,...
        amplitudes_tent,...
        amplitudes_velocities_tent,...
        pressure_amplitudes_tentative,...
        previous_conditions{end}.current_time + dt,...
        dt,...
        new_center_of_mass,...
        new_CM_velocity_times_dt/dt, ...% Divide by dt!
        new_nb_contact_points);

end





function new_probable_next_conditions = advance_conditions_dep(probable_next_conditions, previous_conditions,  ...
    new_nb_contact_points, dt, PROBLEM_CONSTANTS)

    % if previous_conditions <: ProblemConditions; previous_conditions = [previous_conditions]; end
    n = length(previous_conditions); % Determines the order of the method
    if n > 2 || n < 1; throw("Hey!"); end
    extract_symbol = @(jj, field) previous_conditions{jj}.(field);


    nb_harmonics = previous_conditions{end}.nb_harmonics;
    pressure_amplitudes_tentative = probable_next_conditions.pressure_amplitudes;

    % Deformation amplitudes  at all necessary times
    Y = zeros(n+1, nb_harmonics, 2);
    %DA = zeros(n+1, nb_harmonics);
    %DV = zeros(n+1, nb_harmonics);
    % Independent term (pressure term) at all necessary times. 
    % (nb of times used, nb of harmonics, Spatial dimension)
    PB = zeros(1, nb_harmonics, 2); 
    
    if n == 1
        coefs = [-1.0, 1.0];
    elseif n == 2
        rk = dt/previous_conditions{end}.dt;
        ak = (1+2*rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        coefs = [ck, bk, ak]; 
    end
    inverses = PROBLEM_CONSTANTS.ODE_inverse_matrices;
    omegas = PROBLEM_CONSTANTS.omegas_frequencies;
    matrices = PROBLEM_CONSTANTS.ODE_matrices;
    for ii = 2:nb_harmonics %  Do I really need to translate into new coordinates?
        % Translating the problem into new coordinates
        PB(1, ii, :) = inverses(:, :, ii) * ...
            [0; -ii * pressure_amplitudes_tentative(ii)]; 
        %  PB(ii] = -ii * pressure_amplitudes_tentative[ii];
        
        for jj = 1:n
            Y(jj, ii, :) = inverses(:, :, ii) * [ ...
                previous_conditions{jj}.deformation_amplitudes(ii); ...
                previous_conditions{jj}.deformation_velocities(ii)];
            % DA(jj, ii) = previous_conditions(jj).deformation_amplitudes(ii); 
            % DV(jj, ii) = previous_conditions(jj).deformation_velocities(ii);
        end

        %  Extracting tentative solution
        om_i = omegas(ii);
        diagexp = @(dt) expm(diag(-dt .* [-1.0i * om_i, 1.0i * om_i]));
        
        %  This funciton will extract how much time passed 
        time_diference = @(jj) dt + previous_conditions{end}.current_time - previous_conditions{jj}.current_time;

        rhs_vector = zeros(2, 1);
        for jj = 1:n
            rhs_vector = rhs_vector + coefs(jj) * diagexp(time_diference(jj)) * reshape(Y(jj, ii, :), 2, 1);
            % rhs_vector = rhs_vector + coefs[jj] * Y[jj, ii, :]
        end
        Y(end, ii, :) = (dt * reshape(PB(1, ii, :), 2, 1)  - rhs_vector)/ coefs(end);
        Y(end, ii, :) = matrices(:, :, ii) * reshape(Y(end, ii, :), 2, 1);        
    end
    
    amplitudes_tent = Y(end, :, 1);
    amplitudes_velocities_tent = Y(end, :, 2);


    new_CM_velocity_times_dt = - sum(coefs(1:n) .* arrayfun(@(idx) extract_symbol(idx, "center_of_mass_velocity"), 1:n)) * dt;
    
    pressure_amplitudes_tentative = [pressure_amplitudes_tentative, 0];
    Cl = @(l)  l * (l-1) / (2*l-1)     * pressure_amplitudes_tentative(l-1);
    Dl = @(l)  (l+2) * (l+1) / (2*l+3) * pressure_amplitudes_tentative(l+1);
    pressure_amplitudes_tentative = pressure_amplitudes_tentative(1:(end-1));
    
    new_CM_velocity_times_dt = new_CM_velocity_times_dt - dt^2 * PROBLEM_CONSTANTS.froude_nb; 
    %Special case: First harmonics
    new_CM_velocity_times_dt = new_CM_velocity_times_dt - dt^2 * pressure_amplitudes_tentative(1); 
    % General case
    
    for ll = 2:nb_harmonics
        new_CM_velocity_times_dt = new_CM_velocity_times_dt + 3 * dt^2 * ...
            (amplitudes_tent(ll) / (2*ll+1)) * (Cl(ll) - Dl(ll)); 
    end
    %new_CM_velocity_times_dt = new_CM_velocity_times_dt + 3 * dt^2 * ...
    %        (amplitudes_tent(ll) / (2*ll+1)) * (Cl(ll) - 0); 
    

    new_CM_velocity_times_dt = new_CM_velocity_times_dt / coefs(end);
    
    new_center_of_mass = (new_CM_velocity_times_dt - sum(coefs(1:n) .* arrayfun(@(idx) extract_symbol(idx, "center_of_mass"), 1:n)));
    new_center_of_mass = new_center_of_mass /coefs(end);

    new_probable_next_conditions = ProblemConditions( ...
        nb_harmonics,...
        amplitudes_tent,...
        amplitudes_velocities_tent,...
        pressure_amplitudes_tentative,...
        previous_conditions{end}.current_time + dt,...
        dt,...
        new_center_of_mass,...
        new_CM_velocity_times_dt/dt, ...% Divide by dt!
        new_nb_contact_points);

end


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

% function projected_pressure_amplitudes = project_amplitudes(ff, harmonics_qtt, endpoints)
%     projected_pressure_amplitudes = zeros(1, harmonics_qtt);
%     for ii = 1:harmonics_qtt
%         projected_pressure_amplitudes(ii) = integral(@(theta) legendreP(ii, cos(theta)) .* ff(theta) .* sin(theta), endpoints(1), endpoints(2));
%     end
%     % Vectorized vertion
%     %projected_amplitudes = @(ff, harmonics_qtt, endpoints) arrayfun(@(idx) ...
%     %    integral(@(theta) legendreP(ii, cos(theta)) .* ff(theta) .* sin(theta), endpoints(1), endpoints(2)), ...
%     %    1:harmonics_qtt);
% end


%    calculate_exit_angle(::Vector{Float64}, ::Float64)
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


