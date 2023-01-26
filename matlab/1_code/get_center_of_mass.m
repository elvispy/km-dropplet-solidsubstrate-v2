function [CM, CM_velocity] = get_center_of_mass(new_amplitudes, new_velocities, ...
    probable_next_conditions, previous_conditions, PROBLEM_CONSTANTS)
    
    dt = probable_next_conditions.dt;
    n = length(previous_conditions); % Determines the order of the method
    if n > 2 || n < 1; throw("Hey!"); end
    nb_harmonics = PROBLEM_CONSTANTS.nb_harmonics;
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
        CM_velocity = (-dt/PROBLEM_CONSTANTS.froude_nb - sum(coefs(1:n) .* extract_symbol('center_of_mass_velocity')))/coefs(end);
        CM= (dt * CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);    
    else
        switch PROBLEM_CONSTANTS.CM
            case 1 
                % First way: "Exact integration using pressure coefficients
                new_contact_radius = probable_next_conditions.contact_radius; % TODO: MODIFY THIS
                %warning("This version is deprecated!");
                cos_theta = cos(theta_from_cylindrical(new_contact_radius, new_amplitudes)); %theta_from_cylindrical(previous_conditions{end}.contact_radius, new_amplitudes)); 
                %d = 1 + sum(arrayfun(@(idx) (-1)^idx * new_amplitudes(idx), 2:nb_harmonics));
                A = 3/2 * dt * ...
                    sum(probable_next_conditions.pressure_amplitudes .* ...
                        arrayfun(@(idx) integral(@(x) my_legendre(idx, x)./x.^3, -1, cos_theta), 1:nb_harmonics));
                % (Now the B0 = -sum(Bl)  term)
                f = @(x) -0.5 ./(x.^2);
                z = zeta_generator(probable_next_conditions.pressure_amplitudes);
                avg = mean(z(linspace(0, pi/10)));
                A = A + 3/2 * dt * (-avg) * (f(cos_theta) - f(-1));
                %A = - A;

                B = coefs(end)^2/dt;

                C = dt/PROBLEM_CONSTANTS.froude_nb + dot(coefs(1:n), ...
                    (coefs(end)/dt * extract_symbol('center_of_mass') + extract_symbol('center_of_mass_velocity')));

                CM = (-B + sqrt(B^2 - 4*A*C))/(2*A);
                CM_velocity  = (sum(coefs(1:n) .* extract_symbol('center_of_mass')) + coefs(end) * CM)/dt;
                % new_CM_velocity  = (dt * d2zd2t - sum(coefs(1:n) .* extract_symbol('center_of_mass_velocity')))/coefs(end);
                % new_centerofmass = (dt * new_CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);
            case 2 % Second way: Using the small oscillation's hypothesis
                warning("This version is deprecated because the first pressure amplitude is needed!");
                pressure_amplitudes_tentative = [probable_next_conditions.pressure_amplitudes, 0];
                Cl = @(l)  l * (l-1) / (2*l-1)     * pressure_amplitudes_tentative(l-1);
                Dl = @(l)  (l+2) * (l+1) / (2*l+3) * pressure_amplitudes_tentative(l+1);
                pressure_amplitudes_tentative = pressure_amplitudes_tentative(1:(end-1));

                d2zd2t = -1/PROBLEM_CONSTANTS.froude_nb - pressure_amplitudes_tentative(1);%- sum(coefs(1:n) .* extract_symbol("center_of_mass_velocity"), 1:n));

                for hb = 2:nb_harmonics
                    d2zd2t = d2zd2t + 3 * (new_amplitudes(hb) / (2*hb+1)) * (Cl(hb) - Dl(hb));  % THIS IS WRONG IF THE FIRST PRESSURE COEFFICIENT IS DISREGARDED
                end
                CM_velocity  = (dt * d2zd2t - sum(coefs(1:n) .* extract_symbol('center_of_mass_velocity')))/coefs(end);
                CM = (dt * CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);
            case 3
                % Match exactly velocity and center of mass at pi.
                CM_velocity  = sum(arrayfun(@(idx) (-1)^idx * new_velocities(idx), 2:nb_harmonics));
                CM = 1 + sum(arrayfun(@(idx) (-1)^idx * new_amplitudes(idx), 2:nb_harmonics));
            case 4
                % Velocity moves as the amplitudes at pi
                CM_velocity  = sum(arrayfun(@(idx) (-1)^idx * new_velocities(idx), 2:nb_harmonics));
                CM = (dt * CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);
            case 5
                % Zcm accelerates just as the amplitudes at pi
                d2zd2t = - sum(arrayfun(@(idx) (-1)^idx * (idx * probable_next_conditions.pressure_amplitudes(idx) ...
                    + PROBLEM_CONSTANTS.omegas_frequencies(idx)^2 * new_amplitudes(idx)), 2:nb_harmonics));
                CM_velocity  = (dt * d2zd2t - sum(coefs(1:n) .* extract_symbol('center_of_mass_velocity')))/coefs(end);
                CM = (dt * CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);

            case 6
                % Integrate pressure distribution "Exactly on contact area"
                z = zeta_generator(probable_next_conditions.pressure_amplitudes);
                avg = mean(z(linspace(0, pi/10)));
                pb = @(x) (sum(probable_next_conditions.pressure_amplitudes' .* collectPl(probable_next_conditions.nb_harmonics, x), 1) - avg) ./ x.^3;
                endpoints = [-1, cos(theta_from_cylindrical(probable_next_conditions.contact_radius, probable_next_conditions))];
                D = integral(pb , endpoints(1), endpoints(2), 'RelTol', 1e-4);
                A = 3*D*dt/2;
                B = coefs(end)^2/dt;

                C = dt/PROBLEM_CONSTANTS.froude_nb + dot(coefs(1:n), ...
                    (coefs(end)/dt * extract_symbol('center_of_mass') + extract_symbol('center_of_mass_velocity')));
                if A ~= 0
                    CM = (-B + sqrt(B^2 - 4*A*C))/(2*A);
                    CM_velocity  = (sum(coefs(1:n) .* extract_symbol('center_of_mass')) + coefs(end) * CM)/dt;
                else
                    CM_velocity = (-dt/PROBLEM_CONSTANTS.froude_nb - dot(coefs(1:n), extract_symbol('center_of_mass_velocity')))/coefs(end);
                    CM = (dt * CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);
                end
        end % end switch statement
    end
    
    
end