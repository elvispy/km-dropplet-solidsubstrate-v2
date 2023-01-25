function new_probable_next_conditions = pressure_with_flattened_surface(probable_next_conditions, ...
            previous_conditions, PROBLEM_CONSTANTS)
        
    switch PROBLEM_CONSTANTS.PG
        case 1
            if probable_next_conditions.contact_radius < 1e-5
                pressure_guess = probable_next_conditions.pressure_amplitudes/2;
            else

                nb_harmonics = probable_next_conditions.nb_harmonics;

                theta_contact = theta_from_cylindrical(probable_next_conditions.contact_radius, probable_next_conditions);
                M = 5 * nb_harmonics;
                thetas = linspace(theta_contact, pi, M);
                LVAL = collectPl(nb_harmonics, cos(thetas));
                LVAL = LVAL(2:end, :)';
                b = -probable_next_conditions.center_of_mass./cos(thetas) - 1; b = b';

                amplitudes_modified = [0, lsqr(LVAL, b, [], 100, [], [], previous_conditions{end}.deformation_amplitudes(2:end)')'];
            %     g = zeta_generator(probable_next_conditions);
            %     f = @(theta) f_generator(theta, theta_contact, g, probable_next_conditions.center_of_mass);
                %amplitudes_modified = project_amplitudes(f, nb_harmonics, [0, pi, theta_contact], PROBLEM_CONSTANTS, true);
                %amplitudes_modified(1) = 0;        
                pressure_guess = solve_ODE_unkown(amplitudes_modified, nan, probable_next_conditions.dt, previous_conditions, PROBLEM_CONSTANTS);

                % Now we must Calculate B(1) by hand 
            %     n = 2; % TODO UPDATE THIS
            %     dt = probable_next_conditions.dt;
            %     if n == 1
            %         coefs = [-1.0, 1.0];
            %     elseif n == 2
            %         rk = dt/previous_conditions{end}.dt;
            %         ak = (1+2*rk)/(1+rk);
            %         bk = -(1+rk);
            %         ck = rk^2/(1+rk);
            %         coefs = [ck, bk, ak]; 
            %     end
            %     extract_symbol = @(field) arrayfun(@(jj) previous_conditions{jj}.(field), 1:n);
            %    
            %     pressure_amplitudes_tentative = [probable_next_conditions.pressure_amplitudes, 0];
            %     pressure_amplitudes_tentative(1) = 0;
            %     Cl = @(l)  l * (l-1) / (2*l-1)     * pressure_amplitudes_tentative(l-1);
            %     Dl = @(l)  (l+2) * (l+1) / (2*l+3) * pressure_amplitudes_tentative(l+1);
            %     %pressure_amplitudes_tentative = pressure_amplitudes_tentative(1:(end-1));
            % 
            %     d2zd2t = -1/PROBLEM_CONSTANTS.froude_nb;%- sum(coefs(1:n) .* extract_symbol("center_of_mass_velocity"), 1:n));
            %     new_CM_velocity = 0;
            %     for hb = 2:nb_harmonics
            %         d2zd2t = d2zd2t + 3 * (amplitudes_modified(hb) / (2*hb+1)) * (Cl(hb) - Dl(hb)); 
            %         new_CM_velocity  = new_CM_velocity ...
            %             + (-1)^hb * dot(coefs, [arrayfun(@(jj) previous_conditions{jj}.deformation_amplitudes(hb), 1:n), ...
            %             amplitudes_modified(hb)])/dt;
            %     end
                %new_CM_velocity  = (dt * d2zd2t - sum(coefs(1:n) .* extract_symbol('center_of_mass_velocity')))/coefs(end);
                %new_centerofmass = (dt * new_CM_velocity -  sum(coefs(1:n) .* extract_symbol('center_of_mass')))/coefs(end);
            %%case 3

                % JUST TO CHECK THAT B1 is not being used
                pressure_guess(1) = 0;%(d2zd2t - dot(coefs, [extract_symbol('center_of_mass_velocity'), ...
                    %new_CM_velocity])/dt)/(1 - 0.4 * amplitudes_modified(2));
            end
        case 2
            % Let's kill pressure outside
            
    end
    new_probable_next_conditions = probable_next_conditions;
    new_probable_next_conditions.pressure_amplitudes = pressure_guess;
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