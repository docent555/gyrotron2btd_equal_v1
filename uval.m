% Ћинейна€ интерпол€ци€
function u = uval(z, za, r)

persistent i;

n = length(za);

if isempty(i) 
    i = 1;
end
if (i >= n) || (i < 1)
    i = 1;
end

if (z < za(i)) || (z > za(i + 1))
% ! BINARY SEARCH
    i = 1;
    j = n + 1;
    while j > i + 1
        k = fix((i + j)/2);
        if z < za(k)
            j = k;
        end
        if z >= za(k)
            i = k;
        end        
    end
end

u = (r(i)*za(i + 1) - r(i + 1)*za(i))/(za(i + 1) - za(i)) + ...
    (r(i + 1) - r(i))/(za(i + 1) - za(i))*z;
end

