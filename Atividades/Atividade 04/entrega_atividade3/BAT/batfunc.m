function [ x, fval, FOB, Am, rm ] = batfunc( d, N, iter, y, a, MORCEGOS, Linf, Lsup )
% Cálculo da FOB
for m = 1:N
    MORCEGOS(m,d+1) = FOB_ga( MORCEGOS(m,1:d) );
end, clear m
% Atualização do melhor morcego
[index, lin] = find(MORCEGOS(:,end) == min(MORCEGOS(:,end)) );
MELHOR = MORCEGOS(index, :);
% Processo de iteração (parada pr iteração)
t = 1;
r = zeros(1,N);  % Tx de emissão (m,t,val)
A = ones(1,N);   % Amplitude sonora
v = zeros(N, d);   % Velocidae Inicial
while t <= iter
    for m = 1:N
    % Velocidade e deslocamento
    % atualização velocidade
    v(m,:) = v(m,:) + (MELHOR(1:d) - MORCEGOS(m,1:d) )*unifrnd(0,1);
    % Atualização deslocamento
    TEMP(1:d) = MORCEGOS(m,1:d) + v(m,:);
    % Aplicação dos limites
    TEMP(1:d) = limites(TEMP(1:d), Linf , Lsup );
    % Busca local
    if rand < r(m)
        ep = unifrnd(-1, 1, 1, d);
        TEMP(1:d) = MELHOR(1:d) + ep*mean(A);
        % Aplicação dos limites
        TEMP(1:d) = limites(TEMP(1:d), Linf , Lsup);
    end
    % Perturbação em uma dimensão / Avaliação FOB
    if rand < 0.8
        % sorteio dim
        TEMP(1,randsample(1:d,1)) = unifrnd( Linf, Lsup);
    end
    % Aplicação dos limites
    TEMP = limites(TEMP(1:d), Linf ,Lsup );
    % Avaliação FOB
    TEMP(1,d +1) = FOB_ga( TEMP );
    % Busca Global
    if (rand < A(m))||(TEMP(1,end)<=MORCEGOS(m,end))
        MORCEGOS(m,:) = TEMP;  % Atualização do morcego
        r(m) = 1-exp(-y*t);    % atualização da taxa de emissao
        A(m) = a*A(m);         % atualização da amplitude
    end
    % Atualiza MELHOR morcego
    best = find(MORCEGOS(:,end) == min(MORCEGOS(:,end)));
    MELHOR = MORCEGOS(best(1),:);
    FOB(t,:) = MELHOR;
    Am(t,1) = mean(A);
    rm(t,1) = mean(r);
    end
    t = t + 1;
end
x = MELHOR(1,1:d);
fval = MELHOR(1,end);
end

