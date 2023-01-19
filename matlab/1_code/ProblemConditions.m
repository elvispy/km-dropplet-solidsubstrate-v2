function res = ProblemConditions(nb_harmonics, deformation_amplitudes, deformation_velocities, ...
       pressure_amplitudes, current_time, dt, center_of_mass, center_of_mass_velocity, number_contact_points)

res = struct( ...
            "nb_harmonics", nb_harmonics, "deformation_amplitudes", deformation_amplitudes, ...
            "deformation_velocities", deformation_velocities, "pressure_amplitudes", pressure_amplitudes, ...
            "current_time", current_time, "dt", dt, "center_of_mass", center_of_mass, "center_of_mass_velocity", ...
            center_of_mass_velocity, "number_contact_points", number_contact_points);
        
end