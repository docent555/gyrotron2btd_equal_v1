function p = oscill_cmplx(field, ZAxis, Delta, p0)
S1 = griddedInterpolant(ZAxis,field,'spline');
[~, p] = ode45(@(z, p) rhs(z, p, Delta, S1) , ZAxis , p0);
end

