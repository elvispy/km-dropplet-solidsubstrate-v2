function check_fclencurt(f, a, b)
    [x, w] = fclencurt(2^15 + 1, a, b);
    % HYPOTHESIS the weights are updise down
    MI = integral(f, a, b);
    CCI = dot(feval(f, x), w);
    fprintf("RESULTS:\n");
    fprintf("MATLAB INTEGRAL: %.10f\n", MI);
    fprintf("Clenshaw-Curtis Approximation: %.10f\n", CCI);
    
    if abs(MI - CCI) < 1e-7 || abs(MI+CCI) < 1e-7 
        disp("They are pretty close!");
    end
    disp("=======================")
    %assert(norm(integral(f, a, b) + dot(feval(f, x), w)) < 1e-8, sprintf("Didnt work :("));
    %disp("Yes it worked!")
end