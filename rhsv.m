function uv = rhsv(z, pv, delta, a, reidx, imidx)

p = pv(reidx) + 1i*pv(imidx);

u = rhs(z, p, delta, a);

uv = [real(u); imag(u)];

end