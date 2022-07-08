% Solução de FOB com Ant Colony
% Script baseado na aula
% 05/04/2017
% Thainan

clc
clear all
close all 

% Numeros de formigas
Dados_prob.F = 15;
% n° de iterações
Dados_prob.fim_iter = 400;
% taxa de evaporação
Dados_prob.ro = 0.15;
     

%% Dados Circuito
Dados_prob.n_variaveis = 3;
Dados_prob.V1 = 10.6;
Dados_prob.V2 = 10.2;
Dados_prob.V3 = -9.8;
Dados_prob.G1 = 10;
Dados_prob.G2 = 10;
Dados_prob.Is = 110;
Dados_prob.erro = 1e-6;
Dados_prob.I2 = Dados_prob.G2*(Dados_prob.V1 - Dados_prob.V2);
Dados_prob.rm = (Dados_prob.V2 - Dados_prob.V3)/Dados_prob.I2;
Dados_prob.xA = [1 3 5 10 12 15 17 20 22 25 27 30];
Dados_prob.xB = Dados_prob.xA;
Dados_prob.xC = Dados_prob.xA;

%% Coeficientes da FOB
Dados_prob.FOB = [Dados_prob.V2/2,...
                  Dados_prob.V2-Dados_prob.V3,...
                  -Dados_prob.V3/2];

% Restrições da FOB  - Considerando todas as condiçoes ">"
Dados_prob.RESTRICOES = [Dados_prob.V2 0 Dados_prob.V3 Dados_prob.I2];
% Penalidades de cada restrição               
Dados_prob.PENALIDADES = 1e-3;
% inicialização da primeira linha
Solucao.FEROMONIO(1,1:5) = [0 0 0 0 0];
Solucao.FOB       = zeros(1,Dados_prob.F);
%% Formigas iniciais
for k = 1:Dados_prob.F
    Solucao.FORMIGAS(k,:) = [randsample(Dados_prob.xA,1)...
        randsample(Dados_prob.xB,1) ...
        randsample(Dados_prob.xC,1)];
end, clear k

%% Laço principal
Solucao.iter = 0;
h = waitbar(0,'Buscando Solução...');
while Solucao.iter < Dados_prob.fim_iter
%% CALCULO DA FOB   
    Solucao.FOB = Solucao.FORMIGAS*Dados_prob.FOB'...
        - Dados_prob.I2/2; % Aplicação da FOB
    % Restrição 1
    restricao1 = Solucao.FORMIGAS*Dados_prob.RESTRICOES(1, 1:Dados_prob.n_variaveis)';
    teste1 = restricao1 < Dados_prob.RESTRICOES(1,end)*(1 + Dados_prob.erro);
    teste2 = restricao1 > Dados_prob.RESTRICOES(1,end)*(1 - Dados_prob.erro);
    teste = teste1 & teste2;
    teste = teste + 0; % logical2double
    teste(teste==0) = Dados_prob.PENALIDADES;
    Solucao.FOB = Solucao.FOB.*teste;
    clear teste teste1 teste2 teste3 teste4 restricao1 restricao2
%% CALCULO DO FEROMONIO
    Solucao.FOB_1 = Solucao.FOB/sum(Solucao.FOB);
    % cálculo e armazenamento do feromonio
    for k = 1:size(Solucao.FORMIGAS,1)
        % teste se alguma formiga passou lá, e se sim qual posição
        aux = find(Solucao.FEROMONIO(:,1) == Solucao.FORMIGAS(k,1) ...
            & Solucao.FEROMONIO(:,2) == Solucao.FORMIGAS(k,2) ...
            & Solucao.FEROMONIO(:,3) == Solucao.FORMIGAS(k,3));
        % se nenhuma passou teste = 0, se passou teste = posicao;
        if isempty(aux), teste = 0; else teste = aux; end;
        if teste == 0
            % se teste==0, nenhuma formiga passou lá, ou é a primeira iteração
            % se n==0, é a primeira iteração, se n = 1: ultima posição + 1
            %if size(Solucao.FEROMONIO,1) == 1, n = 0;,else, n = 1;, end
            % posição de inserção
            lin = size(Solucao.FEROMONIO,1)+ 1;
            % inserção da posição
            Solucao.FEROMONIO(lin,1:Dados_prob.n_variaveis) = ...
                Solucao.FORMIGAS(k,1:Dados_prob.n_variaveis);
            % inserção do feromonio
            Solucao.FEROMONIO(lin,Dados_prob.n_variaveis+1) = Solucao.FOB_1(k);
            clear n lin
        else
            % se teste diferente de 0, ou seja, alguma formiga ja passou por lá
            % Só atualizar feromonio: antigo mais atual
            
            Solucao.FEROMONIO(teste,Dados_prob.n_variaveis+1) = ...
                Solucao.FEROMONIO(teste,Dados_prob.n_variaveis+1) + Solucao.FOB_1(k);
        end
    end, clear k aux teste
    % normalização
    Solucao.FEROMONIO(:,Dados_prob.n_variaveis+1) = ...
        Solucao.FEROMONIO(:,Dados_prob.n_variaveis+1)/max(Solucao.FEROMONIO(:,Dados_prob.n_variaveis+1));
%% CALCULO DAS PROBABILIDADES
    Solucao.FEROMONIO(isinf(Solucao.FEROMONIO)|isnan(Solucao.FEROMONIO)) = 0;
    Solucao.FEROMONIO(:,end) = Solucao.FEROMONIO(:,Dados_prob.n_variaveis+1)...
        ./sum(Solucao.FEROMONIO(:,Dados_prob.n_variaveis+1));
    Solucao.FEROMONIO(isinf(Solucao.FEROMONIO)|isnan(Solucao.FEROMONIO)) = 0;  
%% ROLETA
    Solucao = roleta3( Dados_prob, Solucao);
%% EVAPORAÇÃO
    Solucao.FEROMONIO(:,Dados_prob.n_variaveis+1) = ...
        (1 - Dados_prob.ro)*Solucao.FEROMONIO(:,Dados_prob.n_variaveis+1);
    % atualização da iteração
    Solucao.iter = Solucao.iter + 1;
    
waitbar(Solucao.iter/Dados_prob.fim_iter,h);
novoTexto = ['Buscando Solução...' num2str(100*Solucao.iter/Dados_prob.fim_iter, '(%0.f)')];
set( get(findobj(h,'type','axes'),'title'), 'string', novoTexto);
end

close(h), clear h

%% Testes da solução
Final.feromonio_max = max(Solucao.FEROMONIO(:,4));
Final.probabilidade_max = max(Solucao.FEROMONIO(:,5));
Final.pos_fer = find( Solucao.FEROMONIO(:,4) == Final.feromonio_max );
Final.pos_prob = find( Solucao.FEROMONIO(:,5) == Final.probabilidade_max );
Final.solucao = Solucao.FEROMONIO(Final.pos_fer, 1:3);
Final.G3 = Final.solucao(1);
Final.G4 = Final.solucao(2);
Final.G5 = Final.solucao(3);
Final.I4 = Solucao.FEROMONIO(Final.pos_fer, 1:3)*Dados_prob.FOB'...
        - Dados_prob.I2/2;
% Teste
Final.A = [ Dados_prob.G1+Dados_prob.G2        -Dados_prob.G2                  0           0
                 -Dados_prob.G2       Dados_prob.G2+Final.G3+Final.G4    -Final.G4      -1
                       0                        -Final.G4             Final.G4+Final.G5  1
         -Dados_prob.rm*Dados_prob.G2 1+Dados_prob.rm*Dados_prob.G2          -1          0];
Final.b = [Dados_prob.Is; 0; 0; 0];
Final.x = [Dados_prob.V1;
     Dados_prob.V2
     Dados_prob.V3
     Final.I4];

Final.B = Final.A*Final.x;

Final.pos_ver = find(Solucao.FEROMONIO(:,1) == 10&...
Solucao.FEROMONIO(:,2) == 30 &...
Solucao.FEROMONIO(:,3) == 10);
Solucao.FEROMONIO(Final.pos_ver,:)
Solucao.FEROMONIO(Final.pos_fer,:)

%% Exibição resultado
clc
fprintf('*******************************\n');
fprintf('              Fim!             \n');
fprintf('*******************************\n');
fprintf ('Resposta:\n\tG3 = %d\n\tG4 = %d\n\tG5 = %d\n\tI4(FOB) = %d\n',...
    Final.solucao, Final.I4 );
fprintf('\tSoluções Visitadas = %d \n', size(Solucao.FEROMONIO,1));
fprintf('\tN° de formigas = %d \n', Dados_prob.F);
fprintf('\tVetor "b"\n');
fprintf('\t\t(%d)\n', int32(Final.B));
fprintf('*******************************\n');



        
        
        


