function [ Solucao ] = roleta3( Dados_prob, Solucao )
% F:        número de formigas      (1,1)
% x:        vetor de possibilidades (1,n)
% PROB:     vetor de probabilidades (1,n)
% formigas: vetor de formigas       (1,F)

% teste da roleta [0 - Vai aleatório (20%)] [1 - Vai na roleta (80%)]
decisao = randsample([0 1], Dados_prob.F, true, [0.2 0.8]);
% Entrada na roleta

for k = 1:Dados_prob.F
    if (decisao(k) == 0 | sum(Solucao.FEROMONIO(:,end)) == 0)% Vai aleatório
        sortear = 1;
        while sortear == 1
            formiga_aleatoria = [randsample(Dados_prob.xA,1) ...
                randsample(Dados_prob.xB,1)...
                randsample(Dados_prob.xC,1)];
            for m = 1:size(Solucao.FEROMONIO,1)
                teste(m,1:Dados_prob.n_variaveis) = Solucao.FEROMONIO(m,1:Dados_prob.n_variaveis) == formiga_aleatoria;
            end
            if sum(teste, 2) == Dados_prob.n_variaveis
                sortear = 1;
            else
                sortear = 0;
            end
        end
            Solucao.FORMIGAS(k,1:Dados_prob.n_variaveis) = formiga_aleatoria;
    else % Vai na roleta
        % sorteio de uma das soluções já existentes no feromonio
        % considerando as probabilidades calculadas
        % aux = (find(Solucao.FEROMONIO(:,end) == 0));
        % Solucao.FEROMONIO(aux,end) = 1e-6;
        solucao = datasample(1:size(Solucao.FEROMONIO,1), 1, 'Weights', Solucao.FEROMONIO(:,end)');
        % atualização da formimga sorteada
        Solucao.FORMIGAS(k,1:Dados_prob.n_variaveis) = Solucao.FEROMONIO(solucao,1:Dados_prob.n_variaveis);
    end
end

end

