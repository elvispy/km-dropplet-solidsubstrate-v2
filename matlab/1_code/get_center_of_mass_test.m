partial_results = zeros(8, 2);
for ii = 1:8
    PROBLEM_CONSTANTS.CM = ii;
    [partial_results(ii, 1), partial_results(ii, 2)] = ...
        get_center_of_mass(new_amplitudes, new_velocities, probable_next_conditions, previous_conditions, PROBLEM_CONSTANTS);
    
end
    