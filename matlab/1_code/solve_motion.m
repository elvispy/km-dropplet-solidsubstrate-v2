% Author: Elvis Agüero
% email: elvisavfc65@gmail.com
% Date: January 13th, 2023.

function solve_motion()

    %solveMotion_2.1
    %Tries to solve the full kinematic match between a dropplet
    % and a solid substrate in vacuum conditions.

    %% Handling default arguments. All units are in cgs.
    
    undisturbed_radius = .1;  % Radius of the undeformed spherical sphere 
    initial_height = Inf;    % Initial position of the sphere center of mass of the sphere (Inf = start barely touching)
    initial_velocity = -.5; % Initial velocity of the sphere 
    initial_amplitudes = Inf; % Initial amplitudes of the dropplet (Default = undisturbed) OBS: First index is A_1
    amplitudes_velocities = [];
    rhoS = 1.0e-3;%%%%            % Sphere's density
    sigmaS = 72.20;%%%%%          % Sphere's Surface Tension
    g = 9.8065e+2;%%%          % Gravitational constant
    harmonics_qtt = 100;      % Number of harmonics to be used 
    nb_pressure_samples = nan;      % Number of intervals in contact radius (NaN = Equal to number of harmonics)
    max_dt = 5e-3;         % maximum allowed temporal time step
    angle_tol = 5/360 * 2 * pi; % Angle tolerance to accept a solution (in radians) 
    spatial_tol = 1e-3;    % Tolerance to accept that dropplet touches the substrate
    simulation_time = 2.0;% Maximum allowed time
    live_plotting = true;     % Whether to plot or not the live results
    


    % Dimensionless Units
    length_unit = undisturbed_radius;
    velocity_unit = abs(initial_velocity);
    time_unit = undisturbed_radius/velocity_unit; % Temporal dimensionless number
    pressure_unit = rhoS * velocity_unit^2;
    froude_nb   = initial_velocity.^2/(g*undisturbed_radius);
    weber_nb    = rhoS * undisturbed_radius * initial_velocity.^2 / sigmaS; % Weber's number of the dropplet
    % reynolds_nb = undisturbed_radius * velocity_unit / nu; % Reynolds' number
    mS = 4*pi*undisturbed_radius^3 * rhoS / 3;
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
    
    initial_velocity = initial_velocity/velocity_unit;

    if length(amplitudes_velocities) ~= harmonics_qtt
        amplitudes_velocities = zeros(1, harmonics_qtt);
    else
        amplitudes_velocities = amplitudes_velocities/velocity_unit;
    end

    tan_tol = tan(angle_tol);

    initial_pressure_coefficients = zeros(1, harmonics_qtt) / pressure_unit; % Just to emphasize the units of these coefficients.

    if max_dt == 0
        dt = 0.01; 
    else 
        dt = max_dt/time_unit; 
    end
    
    % Setting dr
    if isnan(nb_pressure_samples)
        f = @(t) -t^2/(2 * froude_nb) + initial_velocity * t + initial_height;
        spatial_step = sqrt(f(0)^2 - f(dt)^2)/2;
        %spatial_step = 1/harmonics_qtt;
    else
        spatial_step = 1/nb_pressure_samples;
    end
    
    initial_time = 0;
    current_time = initial_time/time_unit;
    final_time = simulation_time/time_unit;
    current_index = 2;%  This integer points to the next available index in variables that are going to 
                      %  export data (index 1 is for initial conditions)
    maximum_index = ceil((final_time - initial_time)/dt) + 4;
    number_of_extra_indexes = 0;

    contact_points = 0;%  Initial number of contact points
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

    f = @(n)  sqrt(n .* (n+2) .* (n-1) / weber_nb);

    omegas_frequencies = f(1:harmonics_qtt)';

    ODE_matrices = zeros(2, 2, harmonics_qtt); % Y' = -PDP^-1 Y + B ==> (exp(tD)*Y)' = e^(tD) P^-1 B;
    ODE_matrices(1, 1, :) =  ones( 1, harmonics_qtt);
    ODE_matrices(1, 2, :) =  ones( 1, harmonics_qtt);
    ODE_matrices(2, 1, :) =  1.0i * omegas_frequencies;
    ODE_matrices(2, 2, :) = -1.0i * omegas_frequencies;

    ODE_inverse_matrices = 1/2 * ones(2, 2, harmonics_qtt);
    ODE_inverse_matrices(1, 2, :) = -0.5i ./ omegas_frequencies;
    ODE_inverse_matrices(2, 2, :) =  0.5i ./ omegas_frequencies;
    
%     syms x;
%     LEGENDRE_POLYNOMIALS = arrayfun( ...
%         @(idx) matlabFunction(legendreP(idx, x)), 1:harmonics_qtt, ...
%         'UniformOutput', false);
%     
%     LEGENDRE_DERIVATIVES = arrayfun( ...
%         @(idx) matlabFunction(diff(legendreP(idx, x), 1)), 1:harmonics_qtt,...
%         'UniformOutput', false);
%     LEGENDRE_DERIVATIVES{1} = @(x) (LEGENDRE_DERIVATIVES{1}());
%     
%     LEGENDRE_SECOND_DERIVATIVES = arrayfun( ...
%         @(idx) matlabFunction(diff(legendreP(idx, x), 2)), 1:harmonics_qtt,...
%         'UniformOutput', false);
%     LEGENDRE_SECOND_DERIVATIVES{1} = @(x) LEGENDRE_SECOND_DERIVATIVES{1}();
%     LEGENDRE_SECOND_DERIVATIVES{2} = @(x) LEGENDRE_SECOND_DERIVATIVES{2}();
%     
%     collectdnPl = @(x) arrayfun(@(idx) LEGENDRE_DERIVATIVES{idx}(x), 1:harmonics_qtt);
%     collectd2nPl = @(x) arrayfun(@(idx) LEGENDRE_SECOND_DERIVATIVES{idx}(x), 1:harmonics_qtt);
%     clear x LEGENDRE_DERIVATIVES LEGENDRE_DERIVATIVES LEGENDRE_SECOND_DERIVATIVES collectdnPl collectd2nPl;
    [nodes, weights] = fclencurt(2^19+1, pi/2, pi);
    
    %names = ["froude_nb" "weber_nb" "omegas_frequencies" "ODE_matrices" "ODE_inverse_matrices"];
    PROBLEM_CONSTANTS = struct("froude_nb", froude_nb, "weber_nb", weber_nb, ...
        "omegas_frequencies", omegas_frequencies, "ODE_matrices", ODE_matrices, ...
        "ODE_inverse_matrices", ODE_inverse_matrices, ...
        "nodes", nodes, "weights", weights, ...
        "pressure_unit", pressure_unit, ...
        "wigner3j", {precomputed_wigner(harmonics_qtt)}, ...
        "DEBUG_FLAG", true);
        
%         "LEGENDRE_POLYNOMIALS", {LEGENDRE_POLYNOMIALS}, ...
%         "LEGENDRE_DERIVATIVES", {LEGENDRE_DERIVATIVES}, ...
%         "collectdnPl", collectdnPl, ...
%         "collectd2nPl", collectd2nPl, ...
%         "pressure_unit", pressure_unit);


    probable_next_conditions = cell(5, 1); % probable_next_conditions = Vector{ProblemConditions}(undef, 5);
    current_conditions = ProblemConditions( ...
        harmonics_qtt, ...
        initial_amplitudes, ...
        amplitudes_velocities, ...
        initial_pressure_coefficients, ...
        current_time, ...
        dt, ...
        initial_height, ...
        initial_velocity, 0);
    
    %current_conditions = ProblemConditions(harmonics_qtt, initial_amplitude, 
    %        amplitudes_velocities, initial_pressure_coefficients, 0.0, dt, 
    %        initial_height, initial_velocity, 0);
    previous_conditions = {current_conditions, current_conditions}; 
    % TODO: Define this array properly to implement BDF2.
    previous_conditions{1}.current_time = previous_conditions{2}.current_time - dt;
    previous_conditions{1}.center_of_mass_velocity = ...
        previous_conditions{2}.center_of_mass_velocity + dt/froude_nb;
    previous_conditions{1}.center_of_mass = ...
        previous_conditions{2}.center_of_mass - previous_conditions{2}.center_of_mass_velocity * dt;
%     currdirr = pwd();
%     if ("julia" in readdir());  cd("julia\\");  end
%     if ("1_code" in readdir()); cd("1_code\\"); end
%     A = match(r"pipeline", file_name)
%     if A !== nothing
%         #A = A[1]; 
%         SAVE_PATH = dirname(file_name);
%     else
%         SAVE_PATH = "../2_pipeline/default/out";
%         file_name = "$SAVE_PATH/$file_name";
%     end
% 
%     if isdir(SAVE_PATH) == false; mkpath(SAVE_PATH); end
%     if (isfile(file_name) == false) && false % TODO: Decide what to export
%         headers_data_frame = DataFrame(
%             ID  = []
%         )
%         CSV.write(file_name, headers_data_frame)
%     end

    % Create logging file


   % % Preparing post-processing
   % TODO: Write post processing variables
   %  Preallocate variables that will be exported (All of them have units!)
   recorded_conditions =cell(maximum_index, 1); % Vector{ProblemConditions}(undef, (maximum_index, )); 
   give_dimensions = @(X) ProblemConditions( ...
       X.nb_harmonics, ...
        X.deformation_amplitudes * length_unit, ...
        X.deformation_velocities * velocity_unit, ...
        X.pressure_amplitudes * (mass_unit * length_unit / (time_unit^2 * length_unit^2)), ...
        X.current_time * time_unit, ...
        X.dt * time_unit, ...
        X.center_of_mass * length_unit, ...
        X.center_of_mass_velocity * velocity_unit, ...
        X.number_contact_points); 

     recorded_conditions{1} = give_dimensions(previous_conditions{end});
%     
%    %  Coefficient of restitution
%     mechanical_energy_in = NaN;
%     mechanical_energy_out = NaN; % TODO: Lab COef of restitution?

    indexes_to_save = zeros(maximum_index, 1); indexes_to_save(1) = 1;
    current_to_save = 2;
    % p = parpool(5);
    while ( current_time < final_time) 
        errortan = Inf * ones(5, 1);
        recalculate = false;

        % First, we try to solve with the same number of contact points
        [probable_next_conditions{3}, errortan(3)] = get_next_step(previous_conditions, contact_points, dt, spatial_step, ...
                spatial_tol, PROBLEM_CONSTANTS);

        if abs(errortan(3)) < tan_tol * 1e-2 % If almost no error, we accept the solution
            current_conditions = probable_next_conditions{3};
            previous_conditions = {previous_conditions{2:end} probable_next_conditions{3}};
        else % If there is some error, we try with diffferent contact points
            % Lets try with one more point
            [probable_next_conditions{4}, errortan(4)] = get_next_step(previous_conditions, contact_points + 1, dt, spatial_step, ...
            spatial_tol, PROBLEM_CONSTANTS);
            % Lets try with one point less
            [probable_next_conditions{2}, errortan(2)] = get_next_step(previous_conditions, contact_points - 1, dt, spatial_step, ...
            spatial_tol, PROBLEM_CONSTANTS);

            if (abs(errortan(3)) > abs(errortan(4)) ||  (abs(errortan(3)) > abs(errortan(2))))
                if abs(errortan(4)) <= abs(errortan(2))
                    % Now lets check with one more point to be sure
                    [~, errortan(5)] = get_next_step(previous_conditions, contact_points + 2, dt, spatial_step, ...
                    spatial_tol, PROBLEM_CONSTANTS);

                    if abs(errortan(4)) < abs(errortan(5)) && errortan(4) < tan_tol
                        % Accept new data 
                        previous_conditions = {previous_conditions{2:end} probable_next_conditions{4}};
                        current_conditions  = probable_next_conditions{4};
                        contact_points      = contact_points + 1;
                    else
                        recalculate = true;
                    end
                else
                    % now lets check if errortan is good enough with one point less
                    [~, errortan(1)] = get_next_step(previous_conditions, contact_points - 2, dt, spatial_step, ...
                    spatial_tol, PROBLEM_CONSTANTS);

                    if abs(errortan(2)) < abs(errortan(1)) && errortan(2) < tan_tol
                        % Accept new data
                        previous_conditions = {previous_conditions{2:end} probable_next_conditions{2}};
                        current_conditions  = probable_next_conditions{2};
                        contact_points      = contact_points - 1;
                    else
                        recalculate = true;
                    end 

                end % End of errortan[4] <= errortan[2]
            else % The same number of contact points may be the best
                if errortan(3) == Inf % ==> All errors are infinity
                    recalculate = true;
                else
                    % Accept new data
                    previous_conditions = {previous_conditions{2:end} probable_next_conditions{3}};
                    current_conditions  = probable_next_conditions{3};
                end
            end
        end % End outer if     while ( current_time < final_time)

        if recalculate == true
                dt = dt/2;
                % Refine time step in index notation 
                iii = iii + 1; jjj = 2 * jjj;
            else
                current_time = current_time + dt; jjj = jjj + 1;
                if mod(jjj, 2) == 0 && grow_dt == true
                    jjj = div(jjj, 2); 
                    iii = iii - 1;
                    % Increase time step
                    dt = 2 * dt;
                    % Decrease the number of time you can make dt bigger
                    grow_dt = false;
                end

                %  TODO: Update Indexes if necessary

                % TODO: % Stored data
                recorded_conditions{current_index} = give_dimensions(current_conditions);
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
                   
                    plot_title = sprintf(" t = %-8.5f, CP = %-8g, \n v_k = %-8.5f, z_k = %-8.5f \n", ...
                       current_time * time_unit, current_conditions.number_contact_points, ...
                            current_conditions.center_of_mass_velocity * velocity_unit, ...
                            current_conditions.center_of_mass* length_unit);
                    plot_condition(1, current_conditions, spatial_step, 10, plot_title);

                else
                    % Do some real-time variable updating here
                end

         end

    end
    
    % Post processing
    indexes_to_save = indexes_to_save(1:(current_to_save-1));
    recorded_conditions = recorded_conditions(indexes_to_save);
    save('simulation.mat', 'recorded_conditions');

 
    
end


        