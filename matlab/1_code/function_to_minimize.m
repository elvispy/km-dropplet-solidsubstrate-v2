function res = function_to_minimize(Xn, previous_conditions, dt, theta_contact, settings)
    % Returns the function to be minimized for the Newton-Raphson Method
    
    Fr = settings.Fr;
    Ps = settings.Ps;
    Flat = settings.Flat;
    weights = settings.weights;
    n = length(previous_conditions); % Determines the order of the method
    if n > 2 || n < 1; throw("Hey!"); end

    N = previous_conditions{end}.nb_harmonics;
    %dt = Xn.dt;
    if n == 1
        coefs = [-1.0, 1.0];
    elseif n == 2
        rk = dt/previous_conditions{end}.dt;
        ak = (1+2*rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        coefs = [ck, bk, ak]; 
    end
    
    extract_symbol = @(jj, field) previous_conditions{jj}.(field);
    previous_velocities = zeros(N, n);
    previous_deformation= zeros(N, n);
    previous_COM        = zeros(1, n);
    previous_COM_vel    = zeros(1, n);
    for ii = 1:n
        previous_velocities(:, ii) = reshape(extract_symbol(ii, 'deformation_velocities'), N, 1);
        previous_deformation(:,ii) = reshape(extract_symbol(ii, 'deformation_amplitudes'), N, 1);
        previous_COM(1, ii)        = extract_symbol(ii, 'center_of_mass');
        previous_COM_vel(1, ii)    = extract_symbol(ii, 'center_of_mass_velocity');
    end
    idxs = 1:(N-1);
    previous_velocities  = previous_velocities(2:end, :);
    previous_deformation = previous_deformation(2:end, :);
    current_deformation  = Xn(idxs);
    %current_deformation  = reshape(current_deformation, N-1, 1);
    current_velocities   = Xn(idxs + (N-1));
    %current_velocities   = reshape(current_velocities, N-1, 1);
    current_pressures    = Xn((2*N-1):(end-2));
    center_of_mass       = Xn(end-1);
    center_of_mass_vel   = Xn(end);
    
    % FIRST ROW BLOCK (Deformation evolution pt1)
    R1 = (sum(coefs .* [previous_deformation, ...
        current_deformation], 2) - dt * current_velocities);

    % Second ROW BLOCK (Deformation evolution pt2 )
    Aidx = (2:N)';
    
    D1N = Aidx .* (Aidx + 2) .* (Aidx - 1);
    
    R2 =  (sum(coefs .* [previous_velocities, current_velocities], 2) ...
         + dt * (Aidx .* current_pressures(3:end) + ...
         D1N .* current_deformation));
     
    if theta_contact < pi
        % THIRD ROW BLOCK (No pressure outside condition)
        theta_i = reshape(linspace(0, theta_contact*0.9, Ps), Ps, 1);
        P1 = collectPl(N, cos(theta_i))';
        P1 = [cos(theta_i), P1];
        R3 = sum(P1 * current_pressures, 2);

        % Fourth ROW BLOCK (flat surface on contact angle condition)
        %L = 10; % This, too, must be unified.
        theta_i2 = reshape(linspace(theta_contact, pi, Flat), Flat, 1); 
        x_i2 = cos(theta_i2);
        P2 = (collectPl(N, x_i2)' + x_i2 .* collectdnPl(N, x_i2)');
        P2 = P2(:, 2:end); % Discard A1
        R4 = sin(theta_i2) .* (1 + sum(P2 * current_deformation, 2));

        % Fifth ROW BLOCK (Center of mass condition)
        oneN = (-1).^(Aidx);
        alternating_sum = sum(oneN .* current_deformation);
        R5 =  (center_of_mass - alternating_sum - 1);

        % Sixth ROW BLOCK (Center of mass diff equation)
        R6 =  sum(coefs .* [previous_COM, center_of_mass] , 2) - dt * center_of_mass_vel; %[zeros(1, 3*N-1), coefs(end), -dt];

        % Seventh ROW BLOCK (center of mass velocity diff equation)
        %oneN = (-1).^(Aidx);
        A = (1 + alternating_sum)^2;
        Cl = manual_intPnxm3(N, -1, cos(theta_contact));
        B = sum(current_pressures .* Cl);

        R7 =  (sum(coefs .* [previous_COM_vel, center_of_mass_vel], 2) - ...
            dt * (-1/Fr + 3/2 * A * B));
            % [-dt * 3 * (current_deformation + oneN) * B, ...
            %zeros(1, N-1), -dt*(3/2 * A) * Cl, 0, coefs(end)];
    else
        R3 = current_pressures;
        R4 = zeros(0, 1);
        R5 = zeros(0, 1);
        R6 = sum(coefs .* [previous_COM, center_of_mass] , 2) - dt * center_of_mass_vel;
        R7 = sum(coefs .* [previous_COM_vel, center_of_mass_vel], 2) + dt/Fr;
    end
    res = weights .* [R1;R2;R3;R4;R5;R6;R7];
    
end