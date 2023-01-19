% Proof of concept that you can approximate 

a = 0;

N = 500;

%bl = zeros(N, 1);

int_kl = zeros(1, N);

kk = 0;
for ll = 1:N
    res = 0;
    for mm = abs(kk-ll):(kk+ll)
       res = res + Wigner3j([kk, ll, mm], [0, 0, 0])^2 * (my_legendre(mm+1, a) ...
           - my_legendre(mm-1, a) - my_legendre(mm+1, -1) + my_legendre(mm-1, -1)); 
    end
    int_kl(kk+1, ll) = res;
end


bl = (((1:N) + 1/2) .* int_kl(1, :))';

escada = @(x) sum(bl .* collectPl(N, x), 1);

b = -1:0.001:1;
c = 1.0 * (b <= a);
d = escada(b) + (a+1)/2;
disp(norm(c-d));
clf;
plot(b, d);
hold on;
scatter(b, c, 'x');