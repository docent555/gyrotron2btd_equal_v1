function s = seval_cmplx(z, n, zax, fre, fim, reb, rec, red, imb, imc, imd)

sre = seval(z, n, zax, fre, reb, rec, red);
sim = seval(z, n, zax, fim, imb, imc, imd);

s = sre + 1i*sim;

end

