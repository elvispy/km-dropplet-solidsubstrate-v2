% Author: Elvis Ag�ero
% email: elvisavfc65@gmail.com
% Date: January 13th, 2023.

function solve_motion_v4()

    %
    %Tries to solve the full kinematic match between a dropplet
    % and a solid substrate in vacuum conditions.
    % THis version uses newton raphson with a coupled set of equations, and
    % the presure vector is 0th indexed

    %% Handling default arguments. All units are in cgs.
    
    undisturbed_radius = .1;  % Radius of the undeformed spherical sphere 
    initial_height = Inf;    % Initial position of the sphere center of mass of the sphere (Inf = start at imminent contact)
    initial_velocity = -10; % Initial velocity of the sphere in cm/s
    initial_amplitudes = Inf; % Initial amplitudes of the dropplet (Default = undisturbed) OBS: First index is A_1
    initial_contact_angle = pi;
    amplitudes_velocities = [];
    rhoS = 0.998;%%%%            % Sphere's density
    sigmaS = 72.20;%%%%%          % Sphere's Surface Tension
    g = 9.8065e+2;%%%          % Gravitational constant
    harmonics_qtt = 50;      % Number of harmonics to be used 
    nb_pressure_samples = nan;      % Number of intervals in contact radius (NaN = Equal to number of harmonics)
    max_dt = 0;         % maximum allowed temporal time step 
    % min_angle = 5/360 * 2 * pi; % Angle tolerance to accept a solution (in radians) 
    spatial_tol = 1e-3;    % Tolerance to accept that dropplet touches the substrate
    angle_tol =  pi/harmonics_qtt;
    simulation_time = 1; % Maximum allowed total time
    live_plotting = true; % Whether to plot or not the live results

    % Dimensionless Units
    length_unit = undisturbed_radius;
    time_unit = sqrt(rhoS * length_unit^3 / sigmaS); %undisturbed_radius/velocity_unit; % Temporal dimensionless number
    velocity_unit = length_unit/time_unit; % abs(initial_velocity);    
    pressure_unit = rhoS * velocity_unit^2;
    froude_nb   = initial_velocity.^2/(g*undisturbed_radius);
    weber_nb    = rhoS * undisturbed_radius * initial_velocity.^2 / sigmaS; % Weber's number of the dropplet
    % reynolds_nb = undisturbed_radius * velocity_unit / nu; % Reynolds' number
    mass_unit = rhoS * length_unit^3;
    
    % % Initial conditions
    % Set dropplet's sphere height initial conditions
    get_initial_height = @(amplitudes) 1 - sum(arrayfun(@(idx) amplitudes(idx) * (-1.0)^(idx), 1:length(amplitudes)));
    
    if length(initial_amplitudes) ~= harmonics_qtt
        initial_amplitudes = zeros(1, harmonics_qtt);
    else
        initial_amplitudes = initial_amplitudes/length_unit;
    end
    if initial_height == Inf
        initial_height = get_initial_height(initial_amplitudes);
    else
        assert(initial_height >= get_initial_height(initial_amplitudes));
        initial_height = initial_height/length_unit;
    end
    
    initial_velocity_adim = initial_velocity/velocity_unit;

    if length(amplitudes_velocities) ~= harmonics_qtt
        amplitudes_velocities = zeros(1, harmonics_qtt);
    else
        amplitudes_velocities = amplitudes_velocities/velocity_unit;
    end

    % tan_tol = tan(min_angle);
    % Zero indexed
    initial_pressure_coefficients = zeros(1, harmonics_qtt+1) / pressure_unit; % Just to emphasize the units of these coefficients.

    if max_dt == 0
        dt = 0.01/ceil(abs(initial_velocity_adim)); 
    else 
        dt = max_dt/time_unit; 
    end
    
    %{
 Setting dr
    if isnan(nb_pressure_samples)
        f2 = @(t) -t^2/(2 * froude_nb) + initial_velocity_adim * t + initial_height;
        spatial_step = sqrt(f2(0)^2 - f2(dt)^2)/2;
        %spatial_step = 1/harmonics_qtt;
    else
        spatial_step = 1/nb_pressure_samples;
    end
    %}
    
    initial_time = 0;
    current_time = initial_time/time_unit;
    final_time = simulation_time/time_unit;
    current_index = 2;%  This integer points to the next available index in variables that are going to 
                      %  export data (index 1 is for initial conditions)
    maximum_index = ceil((final_time - initial_time)/dt) + 4;
    number_of_extra_indexes = 0;

    %contact_points = 0;%  Initial number of contact points
    contact_time = 0;% TODO: Lab contact time and so on?
    grow_dt = false;%  THis variable controls how fast dt can grow
    iii = 0; jjj = 0;%  Indexes to keep track how small is dt compared to max_dt

      
    %%- LEGENDRE_POLYNOMIALS = generate_legen% e_polynomials(harmonics_qtt);
    %%- LPdX = Vector{Function}(undef, harmonics_qtt);
    %%- polynomials_antiderivatives = Matrix{Function}(undef, harmonics_qtt, harmonics_qtt);
    
    
%     @variables x
%     for ii = 1:harmonics_qtt
%         for jj = 1:ii
%             % This data structure at [ii, jj] will return the integral of P_{ii-1} * P_{jj-1}
%             polynomials_antiderivatives[ii, jj] = integrate( LEGENDRE_POLYNOMIALS[ii] * LEGENDRE_POLYNOMIALS[jj])
%         end
%         % This array has the integral of P_{ii-1}/x
%         LPdX[ii] = integrate_poly(LEGENDRE_POLYNOMIALS[ii]);
%     end
    C = 1; % This used to be weber, but now all coeficients are 1.
    f = @(n)  sqrt(n .* (n+2) .* (n-1) / C);
    omegas_frequencies = f(1:harmonics_qtt)';

    PROBLEM_CONSTANTS = struct("froude_nb", froude_nb, "weber_nb", weber_nb, ...
        "nb_harmonics", harmonics_qtt, ...
        "omegas_frequencies", omegas_frequencies, ...
        "spatial_tol", spatial_tol, ...
        "angle_tol", angle_tol, ...
        "pressure_unit", pressure_unit, ...
        "CM", 7, ...
        "PG", 2, ...
        "KILL_OUTSIDE", true, ...
        "DEBUG_FLAG", true); % "wigner3j", {precomputed_wigner(harmonics_qtt)}, ...

    %current_conditions = cell(1, 1); % probable_next_conditions = Vector{ProblemConditions}(undef, 5);
    current_conditions = ProblemConditions_v4( ...
        harmonics_qtt, ...
        initial_amplitudes, ...
        amplitudes_velocities, ...
        initial_pressure_coefficients, ...
        current_time, ...
        dt, ...
        initial_height, ...
        initial_velocity_adim, initial_contact_angle); % Last argument is contact radius
 
    previous_conditions = {current_conditions, current_conditions}; 
    % TODO: Define this array properly to implement BDF2.
    previous_conditions{1}.current_time = previous_conditions{2}.current_time - dt;
    previous_conditions{1}.center_of_mass_velocity = ...
        previous_conditions{2}.center_of_mass_velocity + dt/froude_nb;
    previous_conditions{1}.center_of_mass = ...
        previous_conditions{2}.center_of_mass - previous_conditions{2}.center_of_mass_velocity * dt;
    
    g = @(t, idx) current_conditions.deformation_amplitudes(idx) * cos(f(idx) * t) ...
        + current_conditions.deformation_velocities(idx)/(f(idx)+1e-30) * sin(f(idx) * t); 

    for idx = 1:harmonics_qtt
        previous_conditions{1}.deformation_amplitudes(idx) = g(-dt, idx);
        previous_conditions{1}.deformation_velocities(idx) = (g(0, idx) - g(-2*dt/1000, idx))/(2*dt/1000);
    end

   % % Preparing post-processing
   % TODO: Write post processing variables

   %  Preallocate variables that will be exported (All of them have units!)
   recorded_conditions =cell(maximum_index, 1); % Vector{ProblemConditions}(undef, (maximum_index, )); 
   give_dimensions_v2 = @(X) ProblemConditions_v4( ...
       X.nb_harmonics, ...
        X.deformation_amplitudes * length_unit, ...
        X.deformation_velocities * velocity_unit, ...
        X.pressure_amplitudes * (mass_unit * length_unit / (time_unit^2 * length_unit^2)), ...
        X.current_time * time_unit, ...
        X.dt * time_unit, ...
        X.center_of_mass * length_unit, ...
        X.center_of_mass_velocity * velocity_unit, ...
        X.contact_angle); 

     recorded_conditions{1} = give_dimensions_v2(previous_conditions{end});
%     
%    %  Coefficient of restitution
%     mechanical_energy_in = NaN;
%     mechanical_energy_out = NaN; % TODO: Lab COef of restitution?

    indexes_to_save = zeros(maximum_index, 1); indexes_to_save(1) = 1;
    current_to_save = 2;
    % p = parpool(5);
    while ( current_time < final_time) 
        % First, we try to solve with the same number of contact points

        [current_conditions, errorflag] = get_next_step_v4(previous_conditions, dt, PROBLEM_CONSTANTS);

        if errorflag == true %|| abs(theta_from_cylindrical( current_conditions.contact_radius, current_conditions)...
           %- theta_from_cylindrical(current_conditions.contact_radius, current_conditions)) > 2*PROBLEM_CONSTANTS.angle_tol
            dt = dt/2;
            disp("Se dividio dt por 2");
            % Refine time step in index notation 
            iii = iii + 1; jjj = 2 * jjj;
        else
            previous_conditions = {previous_conditions{2:end} current_conditions};
            
            current_time = current_time + dt; jjj = jjj + 1;
            if mod(jjj, 2) == 0 && grow_dt == true
                jjj = floor(jjj/2); 
                iii = iii - 1;
                % Increase time step
                dt = 2 * dt;
                % Decrease the number of time you can make dt bigger
                grow_dt = false;
            end

            %  TODO: Update Indexes if necessary

            % TODO: % Stored data
            recorded_conditions{current_index} = give_dimensions_v2(current_conditions);
            current_index = current_index + 1; % Point to the next space in memory 

            % If we are in a multiple of max_dt, reset indexes
            if jjj == 2^iii
                jjj = 0;
                grow_dt = true;
                indexes_to_save(current_to_save) = current_index - 1;
                current_to_save = current_to_save + 1;
            else
                number_of_extra_indexes = number_of_extra_indexes + 1;
            end

            if live_plotting == true
                % Do some live plotting here

                plot_title = sprintf(" t = %-8.5f (ms), Contact angle = %-8g deg, \n v_k = %-8.5f cm/s, z_k = %-8.5f cm\n", ...
                   1e+3 * current_time * time_unit, (pi-current_conditions.contact_angle)/pi * 180, ...
                        current_conditions.center_of_mass_velocity * velocity_unit, ...
                        current_conditions.center_of_mass* length_unit);
                plot_condition(1, current_conditions, 3, plot_title);

            else
                % Do some real-time variable updating here
            end

        end

    end
    
    % Post processing
    indexes_to_save = indexes_to_save(1:(current_to_save-1));
    recorded_conditions = recorded_conditions(indexes_to_save);
    save('simulation.mat', 'recorded_conditions', 'length_unit', 'velocity_unit', 'pressure_unit', 'froude_nb', 'weber_nb');

 
    
end