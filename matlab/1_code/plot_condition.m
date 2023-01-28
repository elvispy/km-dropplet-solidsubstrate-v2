function plot_condition(idx, conditions, varargin)
    figure(idx);
    if idx == 2
        set(gcf, 'Position', [780 159 760 586]);
    else
        set(gcf, 'Position', [0   159 760 586]);
    end
    clf;
    hold on;  
    cut = 0.75 * pi;
    sample = [linspace(0, cut, 50), linspace(cut, pi, 50)];
    arrX = sin(sample);
    arrY = cos(sample);
    etas = zeta_generator(conditions);
    if isstruct(conditions)
        height = conditions.center_of_mass;
        plot([-conditions.contact_radius, conditions.contact_radius], [0, 0], 'b', 'LineWidth', 2);
        %yline(, 'y', 'LineWidth', 2);
    else
        height = 1 + sum(arrayfun(@(idx) (-1)^idx * conditions(idx), 2:length(conditions)));
    end
    EtaX = arrayfun(@(angle) sin(angle) * (1+  etas(angle)), sample);
    EtaY = height + arrayfun(@(angle) cos(angle) .* (1+  etas(angle)), sample);
    
    plot( EtaX,EtaY, 'LineWidth',1.5 , 'Color', [.5 .5 .5]);
    plot(-EtaX,EtaY, 'LineWidth',1.5 , 'Color', [.5 .5 .5]);
    
%     gr = -100:100;
%     if nargin > 2
%         dr = varargin{1}; 
%         scatter(dr * gr, 0 * gr, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 2, 'SizeData', 10 * ones(size(gr)));
%     end
    
    scatter(0, 0, 'Marker', 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 2, 'SizeData', 10);
    if isstruct(conditions)
        zps = zeta_generator(conditions.pressure_amplitudes);
        avg = mean(zps(linspace(0, pi/10)));
        ps = @(ang) zps(ang) - avg;
        mps = ps(sample);
        mps(5) = 1;
        quiver(EtaX, EtaY, mps .* (-arrX), mps .* (-arrY));  
    end
    
    if nargin > 2 && isstruct(conditions)
        x = 0.1 + conditions.contact_radius;
    else
        x = 0.05;
    end
    if nargin > 2
        N = varargin{1};
    else
        N = 2;
    end
    xlim([-N*x, N*x]);
    ylim([-2*N/5 * x  , 0.8*2*N*x]);
    %axis equal;
    %ylim([-0.1, 2.5]);
    yline(0, 'r', 'LineWidth', 2);
    
    if nargin > 3 && idx == 2 && isstruct(conditions)
        ss = varargin{2};
        
        s = sprintf("Attempting to fit the solution with contact radius %.3f. \n Previous contact radius: %.3f. Iteration number: %d", ...
            conditions.contact_radius, ss.contact_radius, ss.iteration);
        title(s, 'FontSize', 14);
        legend("Contact Radius");
        %x = xlim;
        y = ylim;
        y = y(2) - (y(2) - y(1))/10;
        x = 0;
        text(x, y, sprintf("v_{cm} = %.3g", conditions.center_of_mass_velocity), 'FontSize', 14);
        %pause(0);
    elseif nargin >= 3
        
        title(varargin{2}, 'FontSize', 14);
        drawnow limitrate;
    end
end