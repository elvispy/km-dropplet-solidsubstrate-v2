function plot_condition(idx, conditions, varargin)
    figure(idx);
    if idx == 2
        set(gcf, 'Position', [780 159 760 586]);
    else
        set(gcf, 'Position', [0   159 760 586]);
    end
    clf; hold on; 
    cut = 0.75 * pi;
    sample = [linspace(0, cut, 50), linspace(cut, pi, 100)];
    arrX = sin(sample);
    arrY = cos(sample);
    etas = zeta_generator(conditions);
    if isstruct(conditions)
        height = conditions.center_of_mass;
    else
        height = 0;
    end
    EtaX = arrayfun(@(angle) sin(angle) * (1+  etas(angle)), sample);
    EtaY = height + arrayfun(@(angle) cos(angle) .* (1+  etas(angle)), sample);
    
    plot( EtaX,EtaY, 'LineWidth',1.5 , 'Color', [.5 .5 .5]);
    plot(-EtaX,EtaY, 'LineWidth',1.5 , 'Color', [.5 .5 .5]);
    
    gr = -100:100;
    if nargin > 2
        dr = varargin{1}; 
        scatter(dr * gr, 0 * gr, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 2, 'SizeData', 10 * ones(size(gr)));
    end
    
    scatter(0, 0, 'Marker', 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 2, 'SizeData', 10);
    if isstruct(conditions)
        zps = zeta_generator(conditions.pressure_amplitudes);
        ps = @(ang) zps(ang) - sum(conditions.pressure_amplitudes);
        quiver(EtaX, EtaY, -ps(sample) .* arrX, -ps(sample) .* arrY);
    end
    
    if nargin > 2 && isstruct(conditions)
        x = dr * max(conditions.number_contact_points, 1);
    else
        x = 0.1;
    end
    if nargin > 3
        N = varargin{2};
    else
        N = 2;
    end
    xlim([-N*x, N*x]);
    ylim([-2*N/5 * x  , 0.8*2*N*x]);
    %axis equal;
    %ylim([-0.1, 2.5]);
    yline(0, 'r', 'LineWidth', 2);
    
    if nargin > 4 && idx == 2 && isstruct(conditions)
        ss = varargin{3};
        
        s = sprintf("Attempting to fit the solution with %d contact points. Previous contact points: %d. \n Iteration Nº: %d", ...
            conditions.number_contact_points, ss.previous_contact_points, ss.iteration);
        title(s, 'FontSize', 14);
        pause(0);
    elseif nargin >= 4
        
        title(varargin{3}, 'FontSize', 14);
        drawnow limitrate;
    end
end