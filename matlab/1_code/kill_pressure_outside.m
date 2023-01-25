function pressure_amplitudes = kill_pressure_outside(probable_next_conditions, theta_max, PROBLEM_CONSTANTS)

    x = cos(theta_max);
    N = probable_next_conditions.nb_harmonics;
    wigner3j = PROBLEM_CONSTANTS.wigner3j;
    
    %Past pressure coeficients 
    Bl = probable_next_conditions.pressure_amplitudes;
    
    % Pressure perturbation coefficients
    Cl = zeros(1, N);
    
    % PRobably all the legendre polynomials evaluation we will need
    Ps = collectPl(2*N+1, x);
    % 
    nps = @(idxs) arrayfun(@(idx) (idx > 0) * Ps(idx - (idx-1) * (idx <= 0)), idxs);
    for k = 1:N
        % First, the B0 term
        Cl(k) = -sum(Bl) * Wigner3j([k 0 k], [0 0 0])^2 * (nps(k-1) - nps(k+1)); % B0 = -sum(Bl) !
        
        % Now the others
        for l = 1:N
            if l > k
                w_m = wigner3j{l, k}.^2;
            else
                w_m = wigner3j{k, l}.^2;
            end
            Cl(k) = Cl(k) + Bl(l) * ( ...
                dot(w_m((abs(k-l)+1):(k+l+1)),-nps((abs(k-l)+1):(k+l+1))) ...
              + dot(w_m((abs(k-l)+1):(k+l+1)), nps((abs(k-l)-1):(k+l-1))) ...
              - w_m(1) * (k==l) ); % When k == l, legendre(-1, 1) == 0, as opposed to all other degrees.
        end
    end
    
    pressure_amplitudes = -(2 * (1:N) + 1/2) .* Cl;
end