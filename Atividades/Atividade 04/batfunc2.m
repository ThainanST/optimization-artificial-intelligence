function [ x, fval ] = batfunc( d, N, iter, y, a, MORCEGOS, Linf, Lsup )

Dados_simu.N = N;
Dados_prob.d = d;
Dados_simu.iter = iter;
Dados_simu.y = y;
Dados_simu.a = a;
Solucao.MORCEGOS = MORCEGOS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passo 3: C�lculo da FOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:Dados_simu.N
    Solucao.MORCEGOS(m,Dados_prob.d+1) = FOB_ga( Solucao.MORCEGOS(m,1:Dados_prob.d) );
end, clear m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passo 4: Atualiza��o do melhor morcego
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[index, lin] = find(Solucao.MORCEGOS(:,end) == min(Solucao.MORCEGOS(:,end)) );
Solucao.MELHOR = Solucao.MORCEGOS(index, :);
clear lin index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passo 5: Processo de itera��o (parada pr itera��o)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 1;
Solucao.r = zeros(1,Dados_simu.N);  % Tx de emiss�o (m,t,val)
Solucao.A = ones(1,Dados_simu.N);   % Amplitude sonora
Solucao.v = zeros(Dados_simu.N, Dados_prob.d);   % Velocidae Inicial

while t <= Dados_simu.iter
    % La�o secund�rio, morcegos
    for m = 1:Dados_simu.N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Passo 7 a 9: Velocidade e deslocamento iniciais
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % atualiza��o velocidade
    Solucao.v(m,:) = Solucao.v(m,:) + ...
        (Solucao.MELHOR(1:Dados_prob.d) - Solucao.MORCEGOS(m,1:Dados_prob.d) )*unifrnd(0,1);
    % Atualiza��o deslocamento
    Solucao.TEMP(1:Dados_prob.d) = Solucao.MORCEGOS(m,1:Dados_prob.d) + Solucao.v(m,:);
    % Aplica��o dos limites
    Solucao.TEMP(1:Dados_prob.d) = limites(Solucao.TEMP(1:Dados_prob.d), Linf , Lsup );
    %% Passo 10 a 12:  Busca local
    if rand < Solucao.r(m)
        % Solucao.TEMP(1:Dados_prob.d) = Solucao.TEMP(1:Dados_prob.d) + 0.001*randn(1,Dados_prob.d);
        ep = unifrnd(-1, 1, 1, Dados_prob.d);
        Solucao.TEMP(1:Dados_prob.d) = Solucao.MELHOR(1:Dados_prob.d) + ep*mean(Solucao.A);
        clear ep
        % Aplica��o dos limites
        Solucao.TEMP(1:Dados_prob.d) = limites(Solucao.TEMP(1:Dados_prob.d), Linf , Lsup);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Passo 13: Perturba��o em uma dimens�o / Avalia��o FOB
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rand < 0.8
    % sorteio dim
    Solucao.TEMP(1,randsample(1:Dados_prob.d,1)) =...
        unifrnd( Linf, Lsup);
    end
    % Aplica��o dos limites
    Solucao.TEMP = limites(Solucao.TEMP(1:Dados_prob.d), Linf ,Lsup );
    % Avalia��o FOB
    Solucao.TEMP(1,Dados_prob.d +1) = FOB_ga( Solucao.TEMP );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Passo 14 a 18: Busca Global
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (rand < Solucao.A(m))||...
            (Solucao.TEMP(1,end)<=Solucao.MORCEGOS(m,end))
        Solucao.MORCEGOS(m,:) = Solucao.TEMP;       % Atualiza��o do morcego
        Solucao.r(m) = 1-exp(-Dados_simu.y*t);      % atualiza��o da taxa de emissao
        Solucao.A(m) = Dados_simu.a*Solucao.A(m);   % atualiza��o da amplitude
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Passo 19: Atualiza MELHOR morcego
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    best = find(Solucao.MORCEGOS(:,end) == min(Solucao.MORCEGOS(:,end)));
    Solucao.MELHOR = Solucao.MORCEGOS(best(1),:);
    clear best
    end, clear m % end for
    t = t + 1;
end, clear t

x = Solucao.MELHOR(1,1:Dados_prob.d);
fval = Solucao.MELHOR(1,end);

end

