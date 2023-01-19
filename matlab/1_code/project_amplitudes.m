function coefs = project_amplitudes(ff, harmonics_qtt, endpoints, PROBLEM_CONSTANTS, flag) % PROBLEM CONSTANTS used to be here 
    %flag = true;
    if flag
        f = @(theta, idx) my_legendre(idx, cos(theta)) .* ff(theta) .* sin(theta);
        try
            f([0 1], 3);
        catch
            f = @(theta, idx) my_legendre(idx, cos(theta))' .* ff(theta) .* sin(theta);
        end
            
        coefs = arrayfun(@(idx) ...
            (2 * idx + 1)/2 * integral(@(theta) f(theta, idx), endpoints(1), endpoints(2), ...
            'RelTol', 5e-3), 1:harmonics_qtt);
    else
        if endpoints(2) > PROBLEM_CONSTANTS.nodes(1) || endpoints(1) < PROBLEM_CONSTANTS.nodes(end)
            [PROBLEM_CONSTANTS.nodes, PROBLEM_CONSTANTS.weights] = fclencurt(2^19+1, endpoints(1), endpoints(2));
            warning("Integration vector recalculated");
        end
        %loc = 1; n = length(PROBLEM_CONSTANTS.nodes);
        
        loc = knnsearch(PROBLEM_CONSTANTS.nodes, endpoints(1));
        disp(loc);
%         while PROBLEM_CONSTANTS.nodes(loc) > endpoints(1) && 2 * loc < n
%             loc = 2 * loc;
%         end
%         if PROBLEM_CONSTANTS.nodes(loc) > endpoints(1)
%             locmax = n;
%             locmin = loc;
%         else
%             locmin = loc/2;
%             locmax = loc;
%         end
% 
%         while locmax - locmin > 1
%             idxmed = round((locmax + locmin)/2);
%             if PROBLEM_CONSTANTS.nodes(idxmed) > endpoints(1)
%                 locmin = idxmed;
%             else
%                 locmax = idxmed;
%             end
%         end
        nodes   = PROBLEM_CONSTANTS.nodes(1:loc);
        weights = PROBLEM_CONSTANTS.weights(1:loc);
        if loc < 1000; warning(srpintf("Too few evaluations for numerical integration: %d", loc)); end
        
%         coefs = zeros(harmonics_qtt, 1)
%         for idx = 1:harmonics_qtt
%             funcs = my_legendre(idx, cos(theta))
%         end
        leg_matrix = collectPl(harmonics_qtt, cos(nodes));
        integral_results =  ff(nodes) .* sin(nodes) .* weights;
        integral_results = leg_matrix .* (integral_results');
        coefs = sum(((1:harmonics_qtt + 1/2))' .* integral_results, 2);
        
%         coefs = arrayfun(@(idx) ...
%             (2*idx+1)/2 * dot(integral(idx, :), weights), 1:harmonics_qtt);
    end
end