function res = JacobianCalculator(Xn, previous_conditions, dt, theta_contact, settings)

    Ps = settings.Ps;
    FLat = settings.Flat;
    weights = settings.weights;
    n = length(previous_conditions); % Determines the order of the method
    if n > 2 || n < 1; throw("Hey!"); end

    N = previous_conditions{end}.nb_harmonics;
    if n == 1
        coefs = [-1.0, 1.0];
    elseif n == 2
        rk = dt/previous_conditions{end}.dt;
        ak = (1+2*rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        coefs = [ck, bk, ak]; 
    end

    %X = [Al, \dot{Al}, Bl, z, \dot{z}]

    % We proceed by newton-raphson

    % FIRST ROW BLOCK (Deformation evolution pt1)
    R1 = [coefs(end) * eye(N-1), -dt * eye(N-1), zeros(N-1, N+1), zeros(N-1, 2)];

    % Second ROW BLOCK (Deformation evolution pt2 )
    Aidx = 2:N;
    %Bidx = 0:N;
    D1N = diag(Aidx .* (Aidx  + 2) .* (Aidx - 1));
    D2N = [zeros(N-1, 2), diag(Aidx)]; % B0 and B1 do not contribute to the motion of the drop
    R2 = [dt*D1N, coefs(end) * eye(N-1), dt * D2N, zeros(N-1, 2)];

    if theta_contact < pi
        % THIRD ROW BLOCK (No pressure outside condition)
        %M = 2;
        theta_i = linspace(0, theta_contact*0.9, Ps)';
        P1 = collectPl(N, cos(theta_i))';
        P1 = [cos(theta_i), P1];
        R3 = [zeros(Ps, 2*(N-1)), P1, zeros(Ps, 2)];

        % Fourth ROW BLOCK (flat surface on contact angle condition)
        %L = 10;
        theta_i2 = linspace(theta_contact, pi, FLat)'; x_i2 = cos(theta_i2);
        P2 = sin(theta_i2) .* (collectPl(N, x_i2)' + x_i2 .* collectdnPl(N, x_i2)');
        P2 = P2(:, 2:end); % Discard A1
        R4 = [P2, zeros(FLat, 2*N+2)];

        % Fifth ROW BLOCK (Center of mass condition)
        oneN = (-1).^(Aidx);
        R5 = [-oneN, zeros(1, 2*N), 1, 0];

        % Sixth ROW BLOCK (Center of mass diff equation)
        R6 = [zeros(1, 3*N-1), coefs(end), -dt];

        % Seventh ROW BLOCK (center of mass velocity diff equation)
        %oneN = (-1).^(Aidx);
        current_deformation = Xn(1:(N-1))';
        current_pressures   = Xn((2*N-1):(end-2))';
        A = (1 + sum(current_deformation .* oneN))^2;
        Cl = manual_intPnxm3(N, -1, cos(theta_contact))';
        B = sum(current_pressures .* Cl);

        R7 = [-dt * 3 * (current_deformation + oneN) * B, ...
            zeros(1, N-1), -dt*(3/2 * A) * Cl, 0, coefs(end)];
    else
        R3 = [zeros(N+1, 2*(N-1)), eye(N+1), zeros(N+1, 2)];
        R4 = zeros(0, 3*N+1);
        R5 = R4;
        R6 = [zeros(1, 3*N-1), coefs(end), -dt];
        R7 = [zeros(1, 3*N-1),   0, coefs(end)];
    end
        
    res = weights .* [R1;R2;R3;R4;R5;R6;R7];
    
end