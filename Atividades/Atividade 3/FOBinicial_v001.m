function [ MORCEGOS ] = FOBinicial_v001( N, d, x1x2y, MORCEGOS )

for k = 1:N
    aux = ([ ones(length(x1x2y),1) x1x2y([1 2],:)']*MORCEGOS(k,1:d)')';
    MORCEGOS(k,d + 1) = mean(abs(aux -  x1x2y(3,:)),2);
end

end

