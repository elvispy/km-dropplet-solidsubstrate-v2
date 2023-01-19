function compare_projection(angles, values, amplitudes)
    % This script compares the projections 
    z = zeta_generator(amplitudes);
    exact_angles = linspace(0, pi, 2000);
    exact_values = arrayfun(@(ang) z(ang), exact_angles);
    
    figure;
    hold on;
    xlim([-0.05, pi+0.05]);
    plot(angles, values, 'r', 'LineWidth', 2)
    plot(exact_angles, exact_values, 'b', 'LineWidth', 2);
    legend('samples', 'with harmonics');

end