function eta = zeta_generator(amplitudes, varargin)
    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    
    order = length(amplitudes);
    
    if nargin == 1
        if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
        eta = @(theta) sum(amplitudes .* collectPl(order, cos(theta)), 1);
    else
        warning("I'm deprecated!");
        LP = varargin{1}.LEGENDRE_POLYNOMIALS;
        eta = @(theta) arrayfun(@(ang) sum(times(amplitudes, arrayfun(@(idx) LP{idx}(cos(ang)), 1:order))), theta);
    end
end