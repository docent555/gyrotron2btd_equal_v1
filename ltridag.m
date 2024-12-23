function u = ltridag(c,a,b,d)
n = length(a);
gam = zeros(n,1);
u = zeros(n,1);
bet = a(n);
if bet == 0.0
    error('tridag_ser: error at code stage 1')
end
u(n) = d(n)/bet;
for j = n - 1:-1:1
    gam(j) = c(j)/bet;
    bet = a(j) - b(j)*gam(j);
    if bet == 0.0
        error('tridag_ser: error at code stage 2')
    end
    u(j) = (d(j) - b(j)*u(j + 1))/bet;
end
for j = 2:n
    u(j) = u(j) - gam(j - 1)*u(j - 1);
end
end

