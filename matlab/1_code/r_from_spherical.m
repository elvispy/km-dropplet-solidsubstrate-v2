function r = r_from_spherical(theta, amplitudes)
    % d = length(amplitudes);
    if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
    r = sin(theta) .* (1 + sum(amplitudes .* collectPl(length(amplitudes), cos(theta)), 1));
%     r = arrayfun(@(angle) ...
%         sin(angle) * (1 + sum(dot(amplitudes, ...
%         arrayfun(@(idx) LP{idx}(cos(angle)), 1:d)))), theta);
end