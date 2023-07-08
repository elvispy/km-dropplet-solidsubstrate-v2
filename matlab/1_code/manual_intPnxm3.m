function val = manual_intPnxm3(N, a, b)
    % Calculates the integral int_a^b P_n(x)/x^3 dx, not manually 
    val = integral(@(x) collectPl(N, x)./x.^3, a, b, 'ArrayValued', true);
    val = [(0.5 /a^2 - 0.5/b^2);val];
end