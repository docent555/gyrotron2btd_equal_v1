function p = oscill_ode(field, Nz, ZAxis, Delta, p0v, reidx, imidx, opts)
% persistent i;
% if isempty(i)
%     i = 1;
% end

% fre = real(field);
% fim = imag(field);
% [reb,rec,red] = cspline(Nz,ZAxis,fre);
% [imb,imc,imd] = cspline(Nz,ZAxis,fim);
% S1 = @(z) seval_cmplx(z, Nz, ZAxis, fre, fim, reb, rec, red, imb, imc, imd);
S1 = griddedInterpolant(ZAxis, field, 'spline');

% if i == 10000
%     for j=1:Nz
%         ss(j) =  S1(ZAxis(j));
%     end
%     figure;
%     plot(ZAxis, imag(ss), ZAxis, imag(field))
%     pause
% end

% opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

[~, pv] = ode89(@(z, p) rhsv(z, p, Delta, S1, reidx, imidx) , ZAxis , p0v, opts);
p = pv(:,reidx) + 1i*pv(:,imidx);

% i = i + 1;
end

