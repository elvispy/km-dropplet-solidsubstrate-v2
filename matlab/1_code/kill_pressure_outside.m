function pressure_amplitudes = kill_pressure_outside(probable_next_conditions, theta_max, PROBLEM_CONSTANTS)

    x = cos(theta_max);
    N = probable_next_conditions.number_contact_points;
    wigner3j = PROBLEM_CONSTANTS.wigner3j;
    
    %Past pressure coeficients 
    Bl = probable_next_conditions.pressure_amplitudes;
    
    % Pressure perturbation coefficients
    Cl = zeros(1, N);
    
    % PRobably all the legendre polynomials evaluation we will need
    Ps = collectPl(2*N+1, x);
    nps = @(idxs) arrayfun(@(idx) (idx > 0) * Ps(idx - (idx-1) * (idx <= 0)), idxs);
    for k = 1:N
        % First, the B0 term
        Cl(k) = -sum(Bl) * Wigner3j([k 0 k], [0 0 0])^2 * (nps(k-1) - nps(k+1));
        
        % Now the others
        for l = 1:N
            wm = wigner3j(k, l).^2;
            Cl(k) = Cl(k) + Bl(l) * ( ...
                dot(wm,-nps((abs(k-l)+1):(k+l+1))) ...
              + dot(wm, nps((abs(k-l)-1):(k+l-1))));
        end
    end
        % Now the special case of k = l
        
    end
    
end