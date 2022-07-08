function [ fob ] = FOB_ga( x )

MORCEGO = [x(1) x(2) x(3)];
x1x2y = [1  1 1 1 1 1 1 1; % A
        -1  0 1 2 4 5 5 6; % B
        -2 -1 0 1 1 2 3 4; % C
        13 11 9 4 11 9 1 -1];
y1 = x1x2y([1 2 3],:)'*MORCEGO';
erro = y1 - x1x2y(end,:)';
fob = sqrt(erro.^2)/length(erro);
fob = mean(fob);

end