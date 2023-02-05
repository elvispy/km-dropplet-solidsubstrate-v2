% Checking delta function distribution


t = 0.8;
b = 0.5;
a = b + (1-b) * t;
x = linspace(b, 1, 2000);
fa = generate_fuction(350, a);
plot(x, fa(x));

function f = generate_fuction(N, a)
    A0 = (1-a)/2;
    a = [1; collectPl(N+1, a)];
    c = (a(1:N) - a(3:end))/2;
    f = @(x) A0 + sum(c.* collectPl(N, x), 1);

end