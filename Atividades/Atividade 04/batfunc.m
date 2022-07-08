function [ x, fval, FOB, Am, rm ] = batfunc( d, N, iter, y, a, MORCEGOS, Linf, Lsup )
% C�lculo da FOB
for m = 1:N
    MORCEGOS(m,d+1) = FOB_ga( MORCEGOS(m,1:d) );
end, clear m
% Atualiza��o do melhor morcego
[index, lin] = find(MORCEGOS(:,end) == min(MORCEGOS(:,end)) );
MELHOR = MORCEGOS(index, :);
% Processo de itera��o (parada pr itera��o)
t = 1;
r = zeros(1,N);  % Tx de emiss�o (m,t,val)
A = ones(1,N);   % Amplitude sonora
v = zeros(N, d);   % Velocidae Inicial
while t <= iter
    for m = 1:N
    % Velocidade e deslocamento
    % atualiza��o velocidade
    v(m,:) = v(m,:) + (MELHOR(1:d) - MORCEGOS(m,1:d) )*unifrnd(0,1);
    % Atualiza��o deslocamento
    TEMP(1:d) = MORCEGOS(m,1:d) + v(m,:);
    % Aplica��o dos limites
    TEMP(1:d) = limites(TEMP(1:d), Linf , Lsup );
    % Busca local
    if rand < r(m)
        ep = unifrnd(-1, 1, 1, d);
        TEMP(1:d) = MELHOR(1:d) + ep*mean(A);
        % Aplica��o dos limites
        TEMP(1:d) = limites(TEMP(1:d), Linf , Lsup);
    end
    % Perturba��o em uma dimens�o / Avalia��o FOB
    if rand < 0.8
        % sorteio dim
        TEMP(1,randsample(1:d,1)) = unifrnd( Linf, Lsup);
    end
    % Aplica��o dos limites
    TEMP = limites(TEMP(1:d), Linf ,Lsup );
    % Avalia��o FOB
    TEMP(1,d +1) = FOB_ga( TEMP );
    % Busca Global
    if (rand < A(m))||(TEMP(1,end)<=MORCEGOS(m,end))
        MORCEGOS(m,:) = TEMP;  % Atualiza��o do morcego
        r(m) = 1-exp(-y*t);    % atualiza��o da taxa de emissao
        A(m) = a*A(m);         % atualiza��o da amplitude
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

