function [ fob ] = FOB_v001( x1x2y, d, BEST )


aux = ([ ones(length(x1x2y),1) x1x2y([1 2],:)']*BEST(1,1:d)')';
fob = mean(abs(aux - x1x2y(3,:)),2);



end