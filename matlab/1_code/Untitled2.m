theta = linspace(0, pi, 1000);
p1 = zeta_generator(new_probable_next_conditions.pressure_amplitudes);
p2 = zeta_generator(npw);
figure(6)
plot(theta, p1(theta));
hold on;
plot(theta, p2(theta));
x= solve_ODE_unkown(nan, npw, probable_next_conditions.dt, previous_conditions, PROBLEM_CONSTANTS);
y = solve_ODE_unkown(nan, new_probable_next_conditions.pressure_amplitudes, probable_next_conditions.dt, previous_conditions, PROBLEM_CONSTANTS);
size(x)