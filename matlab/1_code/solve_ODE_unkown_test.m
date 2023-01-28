N = 20;
IDX = 15;
M = 200;
dt = 1e-5;
initial_deformation = zeros(1, N);
initial_velocities  = ones(1, N);

f = @(n)  sqrt(n .* (n+2) .* (n-1) / 1);
omegas_frequencies = f(1:N)';
PROBLEM_CONSTANTS = struct("omegas_frequencies", omegas_frequencies);

g = @(t, idx) initial_deformation(idx) * cos(f(idx) * t) + initial_velocities(idx)/(f(idx)+1e-30) * sin(f(idx) * t); 

initial_condition = struct("dt", dt, "nb_harmonics", N, ...
    "deformation_velocities", initial_velocities, "deformation_amplitudes", initial_deformation);
initial_condition_2 = initial_condition;
for idx = 1:N
    initial_condition_2.deformation_amplitudes(idx) = g(-dt, idx);
    initial_condition_2.deformation_velocities(idx) = (g(0, idx) - g(-2*dt/1000, idx))/(2*dt/1000);
end

previous_conditions = {initial_condition_2 initial_condition};
conditions = cell(M, 1); conditions{1} = initial_condition;


for ii = 2:M
    [def, vel] = solve_ODE_unkown(nan, zeros(1, N), dt, previous_conditions, PROBLEM_CONSTANTS);
   
    new_conditions = struct("dt", dt, "nb_harmonics", N, ...
     "deformation_velocities", vel, "deformation_amplitudes", def);

    previous_conditions = {previous_conditions{end} new_conditions};
    conditions{ii} = new_conditions; 
end


t = linspace(0, dt*M, M);

defs = zeros(1, M);
real = zeros(1, M);
for ii = 1:M
    defs(ii) = conditions{ii}.deformation_amplitudes(IDX);
    real(ii) = g(t(ii), IDX);
end

clf; figure(3);
plot(t, defs, 'r');
hold on;
plot(t, real, 'b');
legend(["Numerical" "Analytic"]);
title(sprintf("Amplitude %d", IDX));


%% Test that solve_ODE is an idempotent operator.

A = rand(1, N);
[B, C] = solve_ODE_unkown(nan, A, dt, previous_conditions, PROBLEM_CONSTANTS);
D = solve_ODE_unkown(B, nan, dt, previous_conditions, PROBLEM_CONSTANTS);

E = rand(1, N);
F = solve_ODE_unkown(E, nan, dt, previous_conditions, PROBLEM_CONSTANTS);
G = solve_ODE_unkown(nan, F, dt, previous_conditions, PROBLEM_CONSTANTS);

fprintf("La norma 2 de las presiones     es %.3g \n", norm(A-D));

fprintf("La norma 2 de las deformaciones es %.3g \n", norm(E-G));
