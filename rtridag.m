function u = rtridag(c, a, b, d)
n = length(a);
gam = zeros(n,1);
u = zeros(n,1);
bet = a(1);
if bet == 0.0
    error('tridag: error at code stage 1');
end
u(1) = d(1)/bet;
for j = 2:n
    gam(j) = b(j - 1)/bet;
    bet = a(j) - c(j - 1)*gam(j);
    if bet == 0.0
        error('tridag: error atu code stage 2');
    end
    u(j) = (d(j) - c(j - 1)*u(j - 1))/bet;
end
for j = n - 1:-1:1
    u(j) = u(j) - gam(j + 1)*u(j + 1);
end
end

