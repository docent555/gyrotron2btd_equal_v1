function s = seval(u, n, x, y, b, c, d)
%    x(n), y(n), b(n), c(n), d(n)
% !
% !THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION
% !
% !SEVAL = Y(I)+B(I)*(U-X(I)) + C(I)*(U-X(I)))**2 + D(I)*(U-X(I))**3
% !
% !WHERE X(I) .LT. U .LT. X(I + 1). USING HORNER'S RULE
% !
% !IF U .LT. X(1), THEN THE VALUE 1 = 1 IS TAKEN.
% !IF U .GE. X(N), THEN THE VALUE I = N IS TAKEN.
% !
% !INPUT..
% !
% !N = THE NUMBER OF DATA POINTS
% !U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
% !X, Y = THE ARRAYS OF DATA ABSCISSAS AND ORD1NATES
% !B, C, D = ARRAYS OF SPLINE COEFFICIENTS, COMPUTED BY SPLINE SUBROUTINE
% !
% !IF U IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
% !BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
% !
persistent i;

if isempty(i) 
    i = 1;
end
if (i >= n) || (i < 1)
    i = 1;
end

if (u < x(i)) || (u > x(i + 1))
% ! BINARY SEARCH
    i = 1;
    j = n + 1;
    while j > i + 1
        k = fix((i + j)/2);
        if u < x(k)
            j = k;
        end
        if u >= x(k)
            i = k;
        end        
    end
end
% !
% ! EVALUATE SPLINE
% !
dx = u - x(i);
s = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)));
end
