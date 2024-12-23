function u = rhs(z, p, delta, a)

u = -a(z) - 1i*p.*(delta - 1.0D0 + abs(p).^2);

end
