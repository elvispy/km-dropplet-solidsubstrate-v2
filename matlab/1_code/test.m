% Checking delta function distribution


t = 0.8;
b = 0.5;
a = b + (1-b) * t;
x = linspace(b, 1, 2000);
fa = generate_fuction(350, a);
%plot(x, fa(x));



%% Testing function which determines contat radius

angle_tol = pi/1000;
N = 500;
d = 0.1;
A0 = (1-cos(d))/2;
a = [1; collectPl(N+1, cos(d))];
c = (a(1:N) - a(3:end))/2;
c = c .* (-mod(1:N, 2) + mod((1:N) + 1, 2))';
pc = struct("pressure_amplitudes", c);
z = zeta_generator(pc.pressure_amplitudes);
z = @(ang) z(ang) + A0;
xx = linspace(0, pi, 2000);
plot(xx, z(xx));
xang = pressure_angle(pc, nan, angle_tol);

function c_angle = pressure_angle(probable_next_conditions, ~, angle_tol)
%%     pressure_amplitudes_tentative = [probable_next_conditions.pressure_amplitudes, 0];
%     Cl = @(l)  l * (l-1) / (2*l-1)     * pressure_amplitudes_tentative(l-1);
%     Dl = @(l)  (l+2) * (l+1) / (2*l+3) * pressure_amplitudes_tentative(l+1);
%     pressure_amplitudes_tentative = pressure_amplitudes_tentative(1:(end-1));
% 
%     inst_force = - pressure_amplitudes_tentative(1);%- sum(coefs(1:n) .* extract_symbol("center_of_mass_velocity"), 1:n));
% 
%     for hb = 2:nb_harmonics
%         inst_force = inst_force + 3 * (new_amplitudes(hb) / (2*hb+1)) * (Cl(hb) - Dl(hb));  % THIS IS WRONG IF THE FIRST PRESSURE COEFFICIENT IS DISREGARDED
%     end
%     
%     inst_force = 4*pi/3 * inst_force; % We take out the dimensionless mass
    %% Code
    z = zeta_generator(probable_next_conditions.pressure_amplitudes);
    avg = mean(z(linspace(0, pi/10)));
    z = @(ang) z(ang) - avg;
    inst_force = integral(z, pi/2, pi, 'RelTol',1e-5);

    c_angle = pi;
    d_angle = pi/10;
    F = 0;
    while d_angle > angle_tol    
        c_angle = c_angle - d_angle;
        Fold = F;
        F = F + integral(@(ang) z(ang), c_angle, c_angle + d_angle, 'RelTol',1e-3);
        if F > 0.99 * inst_force && d_angle < angle_tol
            break;
        elseif F > 0.99 * inst_force
            c_angle = c_angle + d_angle;
            d_angle = d_angle/2;
            F = Fold;
        end
    end

end



function f = generate_fuction(N, a)
    A0 = (1-a)/2;
    a = [1; collectPl(N+1, a)];
    c = (a(1:N) - a(3:end))/2;
    f = @(x) A0 + sum(c.* collectPl(N, x), 1);

end