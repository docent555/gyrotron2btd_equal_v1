function [b, c, d] = spline(n, x, y) %#codegen
% x = zeros(n,1);
% y = zeros(n,1);
b = zeros(n,1);
c = zeros(n,1);
d = zeros(n,1);
% !  THE COEFFICIENTS B(I), C(I) AND D(I) ARE CALCULATED, 1=1,
% !  2, ..., N, FOR CUBIC INTERPOLATION SPLINE
% !
% !  S(X) = Y(I)+B(I)*(X-X(I)) + C(I)*(X-X(I))**2 +
% !  -fD(I)*(X - X(I))**3
% !
% !  FOR X(I) .LE. X .LE. X(I+1)
% !
% !  INPUT INFORMATION..
% !
% !  N = NUMBER OF SPECIFIED POINTS OR NODES (N .GE. 2)
% !  X = ABSCISSUE OF NODES IN STRICTLY INCREASING ORDER
% !  Y = ORDINATES OF NODES
% !
% !  OUTPUT...
% !
% !  B, C, D = ARRAYS OF SPLINE COEFFICIENTS DEFINITED ABOVE.
% !
% !  IF YOU DESIGNATE THE DIFFERENTIATION SYMBOL BY P, THEN
% !
% !  Y(I)= S(X(I))
% !  B(I) = SP(X(I))
% !  C(I) = SPP(X(I))/2
% !  D(I) = SPPP(X(I))/6 (RIGHT HAND DERIVATIVE)
% !
% !  USING THE ACCOMPANYING SEVAL FUNCTION SUBROUTINE
% !  YOU CAN CALCULATE SPLINE VALUES.
% !
nm1 = n - 1;
if (n < 2)
    return
end
if (n < 3)
    b(1) = (y(2) - y(1))/(x(2) - x(1));
    c(1) = 0.;
    d(1) = 0.;
    b(2) = b(1);
    c(2) = 0.;
    d(2) = 0.;
    return
end
% !
% ! BUILD A TRIDIAGONAL SYSTEM
% ! B = DIAGONAL, O = OVERDIAGONAL, C = RIGHT PARTS.
% !
d(1) = x(2) - x(1);
c(2) = (y(2) - y(1))/d(1);
for i = 2:nm1
    d(i) = x(i + 1) - x(i);
    b(i) = 2.*(d(i - 1) + d(i));
    c(i + 1) = (y(i + 1) - y(i))/d(i);
    c(i) = c(i + 1) - c(i);
end
% !
% ! BOUNDARY CONDITIONS. THIRD DERIVATIVES AT POINTS
% ! X(1) AND X(N) ARE CALCULATED USING DIVISIONED
% ! DIFFERENCES
% !
b(1) = -d(1);
b(n) = -d(n - 1);
c(1) = 0.;
c(n) = 0.;
if (n > 3)
    c(1) = c(3)/(x(4) - x(2)) - c(2)/(x(3) - x(1));
    c(n) = c(n - 1)/(x(n) - x(n - 2)) - c(n - 2)/(x(n - 1) - x(n - 3));
    c(1) = c(1)*d(1).^2/(x(4) - x(1));
    c(n) = -c(n)*d(n - 1).^2/(x(n) - x(n - 3));
end
% !
% ! STRAIGHT RUN
% !
for i = 2:n
    t = d(i - 1)/b(i - 1);
    b(i) = b(i) - t*d(i - 1);
    c(i) = c(i) - t*c(i - 1);
end
% !
% ! REVERSE SUBSTITUSTION
% !
c(n) = c(n)/b(n);
for ib = 1:nm1
    i = n - ib;
    c(i) = (c(i) - d(i)*c(i + 1))/b(i);
end
% !
% ! C(I) NOW STORES THE VALUE OF SIGMA(I), DEFINED
% ! IN #4.4.
% !
% ! CALCULATE COEFFICIENTS OF POLYNOMIALS
% !
b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n));
for i = 1:nm1
    b(i) = (y(i + 1) - y(i))/d(i) - d(i)*(c(i + 1) + 2.*c(i));
    d(i) = (c(i + 1) - c(i))/d(i);
    c(i) = 3.*c(i);
end
c(n) = 3.*c(n);
d(n) = d(n - 1);
return
end

