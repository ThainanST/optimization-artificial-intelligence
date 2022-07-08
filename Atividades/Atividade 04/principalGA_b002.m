% 02/05/2017
% Técnicas Inteligentes - Bat Algorithm vc Algoritmo Genético
% Prof. Ivo Chaves

clc
clear all
close all

teste = [
10	 50 	30	1e-9	0.05	0.70
40	 50 	30	1e-9	0.05	0.70
70	 50 	30	1e-9	0.05	0.70
100	 50 	30	1e-9	0.05	0.70
40	30	30	1e-9	0.05	0.70
40	50	30	1e-9	0.05	0.70
40	75	30	1e-9	0.05	0.70
40	 50 	10	1e-9	0.05	0.70
40	 50 	25	1e-9	0.05	0.70
40	 50 	50	1e-9	0.05	0.70
40	 50 	30	1e-3	0.05	0.70
40	 50 	30	1e-6	0.05	0.70
40	 50 	30	1e-9	0.05	0.70
40	 50 	30	1e-9	0.0	0.70
40	 50 	30	1e-9	0.05	0.70
40	 50 	30	1e-9	0.10	0.70
40	 50 	30	1e-9	0.05	0.50
40	 50 	30	1e-9	0.05	0.70
40	 50 	30	1e-9	0.05	0.90];

time = 0;
Nsimu = 50;  % N° de simulações
solInit = 0; % 1 = criar; 0 carregar
d = 3;       % N° de dimensões
Linf = 3;  % Limite inferior 3;%
Lsup = 10;   % Limite superior
SOLUCOES0 = unifrnd(Linf,Lsup, 150, d);

h = waitbar(0,'Buscando Soluções... (0)');
for t = 1:size(teste,1)
    % Testes
    N   = teste(t,1);           % N° de individuos
    ger = teste(t,2);           % N° de iterações
    MSG = teste(t,3);           % MaxStallGenerations
    FT  = teste(t,4);           % FunctionTolerance
    EL  = ceil(teste(t,5)*N);   % EliteCount
    CO  = teste(t,6);           % CrossoverFraction   
    
    %% Opções Simulação
    % https://www.mathworks.com/help/gads/genetic-algorithm-options.html
    % Plot Options ('PlotFcn')
    % Best fitness (@gaplotbestf)
    % Best individual (@gaplotbestindiv)
    % Average Pareto distance (@gaplotparetodistance)
    % Population Options
    % Population size (PopulationSize)
    % Initial population (InitialPopulationMatrix)
    % Fitness Scaling Options ('FitnessScalingFcn')
    % Rank (@fitscalingrank)
    % Proportional (@fitscalingprop)
    % Top (@fitscalingtop)
    % Shift linear (@fitscalingshiftlinear)
    % Selection Options ('SelectionFcn')
    % Stochastic uniform (@selectionstochunif)
    % Remainder (@selectionremainder)
    % Uniform (@selectionuniform)
    % Roulette (@selectionroulette)
    % Tournament (@selectiontournament) // options = optimoptions('ga','SelectionFcn', {@selectiontournament,size})
    % Reproduction Options
    % Elite count (EliteCount)
    % Mutation Options ('MutationFcn')
    % {@mutationgaussian, scale, shrink}
    % {@mutationuniform, rate}
    % Crossover Options ('CrossoverFcn')
    % Scattered (@crossoverscattered)
    % point (@crossoversinglepoint)
    % Two point (@crossovertwopoint)
    % Intermediate ({@crossoverintermediate, ratio})
    % Heuristic {@crossoverheuristic,ratio}
    % Arithmetic (@crossoverarithmetic)
    % Stopping Criteria Options
    % Generations (MaxGenerations)
    % Time limit (MaxTime)
    % Fitness limit (FitnessLimit)
    % Stall generations (MaxStallGenerations)
    % Stall time limit (MaxStallTime)
    % Function tolerance (FunctionTolerance)
    % Constraint tolerance (ConstraintTolerance)
    % Output Function Options
    % Display to Command Window Options
    
    %% Estrutura de opções
    SOLUCOES = SOLUCOES0(1:N, :);
    opt = optimoptions(@ga, ... 'PlotFcn', {@gaplotbestf @gaplotbestindiv}, ...
        'InitialPopulationMatrix', SOLUCOES,...
        'PopulationSize', N, ... (%%%)
        'MaxGenerations', ger, ... (%%%)
        'MaxStallGenerations', MSG, ... (%%%)
        'FunctionTolerance', FT, ... (%%%)
        'EliteCount', EL, ... (%%%)
        'CrossoverFraction', CO,... (%%%)
        'SelectionFcn', @selectionroulette,...
        'CrossoverFcn', @crossoversinglepoint, ...
        'MutationFcn', {@mutationadaptfeasible, 1/100}); % mutationuniform mutationadaptfeasible
    % optimoptions gaoptimset     'Display', 'off',...
    
    %% Simulações
    for k = 1:Nsimu
        tic
        [x,fval, flag, output] = ga(@FOB_ga,d,[],[],[],[],Linf*[1 1 1],Lsup*[1 1 1],[],opt);
        AG(k,1:d) = x;
        AG(k,d+1) = fval;
        AG(k,d+2) = toc;
        AG(k,d+3) = flag;
    end
    
    [a, b] = sort(AG(:,d+1));
    AG_sort_fval = AG(b,:);
    
    [a, b] = sort(AG(:,end));
    AG_sort_time = AG(b,:);
    
    Resul = [AG_sort_fval([1 end],:);
        AG_sort_time([1 end],:)];
    Val_medio = mean([AG(:,4) AG(:,5)]);
    Resultado.val{t,1} = [Resul ones(4,1)*Val_medio(1) ones(4,1)*Val_medio(2)];
    
    waitbar(t/length(teste),h);
    ntext = ['Buscando Solução...' num2str(100*t/length(teste), '(%0.f)')];
    set( get(findobj(h,'type','axes'),'title'), 'string', ntext);
    clear ntext
    clc
end
close(h)

for k = 1:size(teste)
    Resultado.final(k,:) = Resultado.val{k}(1,:);
end


%% Testes 1,4
figure
aux = 1:4;
yyaxis left
plot(aux, Resultado.final(aux,4), 'b' ), hold on, grid on
plot(aux, Resultado.final(aux,7), '--k' )
ylabel('FOB');
xlabel('Teste');
yyaxis right
plot(aux, Resultado.final(aux,8));
ylabel('Tempo de simulação, s');
title('Testes 1 a 4: Variação de indivíduos, [3, 10]');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')

% Evolução da fob
SOLUCOES0 = unifrnd(3,10, 150, 3);
% Testes
N   = 100        % N° de individuos
ger = 50;            % N° de iterações
MSG = 30;             % MaxStallGenerations
FT  = 1e-9;           % FunctionTolerance
EL  = ceil(0.05*N);   % EliteCount
CO  = 0.70;            % CrossoverFraction
Linf = 3;
Lsup =10;
% Estrutura de opções
SOLUCOES = SOLUCOES0(1:N, :);
opt = optimoptions(@ga, ...
    'PlotFcn', {@gaplotbestf},... @gaplotbestindiv}, ...
    'InitialPopulationMatrix', SOLUCOES,...
    'PopulationSize', N, ... (%%%)
    'MaxGenerations', ger, ... (%%%)
    'MaxStallGenerations', MSG, ... (%%%)
    'FunctionTolerance', FT, ... (%%%)
    'EliteCount', EL, ... (%%%)
    'CrossoverFraction', CO,... (%%%)
    'SelectionFcn', @selectionroulette,...
    'CrossoverFcn', @crossoversinglepoint, ...
    'MutationFcn', {@mutationadaptfeasible, 1/100}); % mutationuniform mutationadaptfeasible
% optimoptions gaoptimset     'Display', 'off',...
[x,fval, flag] = ga(@FOB_ga,d,[],[],[],[],Linf*[1 1 1],Lsup*[1 1 1],[],opt);


%% variaçãõ de gerações
figure
aux = 5:7;
yyaxis left
plot(aux, Resultado.final(aux,4), 'b' ), hold on, grid on
plot(aux, Resultado.final(aux,7), '--r' )
title('Testes 5 a 7: Variação gerações, [3, 10]');
xlabel('Teste');
ylabel('FOB');
yyaxis right
plot(aux, Resultado.final(aux,8));
ylabel('Tempo de simulação, s');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')

SOLUCOES0 = unifrnd(3,10, 150, 3);
% Testes
N   = 40;        % N° de individuos
ger = 30;            % N° de iterações
MSG = 30;             % MaxStallGenerations
FT  = 1e-9;           % FunctionTolerance
EL  = ceil(0.05*N);   % EliteCount
CO  = 0.70;            % CrossoverFraction
Linf = 3;
Lsup =10;
% Estrutura de opções
SOLUCOES = SOLUCOES0(1:N, :);
opt = optimoptions(@ga, ...
    'PlotFcn', {@gaplotbestf},... @gaplotbestindiv}, ...
    'InitialPopulationMatrix', SOLUCOES,...
    'PopulationSize', N, ... (%%%)
    'MaxGenerations', ger, ... (%%%)
    'MaxStallGenerations', MSG, ... (%%%)
    'FunctionTolerance', FT, ... (%%%)
    'EliteCount', EL, ... (%%%)
    'CrossoverFraction', CO,... (%%%)
    'SelectionFcn', @selectionroulette,...
    'CrossoverFcn', @crossoversinglepoint, ...
    'MutationFcn', {@mutationadaptfeasible, 1/100}); % mutationuniform mutationadaptfeasible
% optimoptions gaoptimset     'Display', 'off',...
[x,fval, flag] = ga(@FOB_ga,d,[],[],[],[],Linf*[1 1 1],Lsup*[1 1 1],[],opt);

%% variaçãõ de gerações estagnadas
figure
aux = 8:10;
yyaxis left
plot(aux, Resultado.final(aux,4), 'b' ), hold on, grid on
plot(aux, Resultado.final(aux,7), '--r' )
title('Testes 8 a 10: Variação gerações estagnadas, [3, 10]');
xlabel('Teste');
ylabel('FOB');
yyaxis right
plot(aux, Resultado.final(aux,8));
ylabel('Tempo de simulação, s');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')

figure
aux = 11:13;
yyaxis left
plot(aux, Resultado.final(aux,4), 'b' ), hold on, grid on
plot(aux, Resultado.final(aux,7), '--r' )
title('Testes 11 a 13: Variação da tolerância estagnação, [-10, 10]');
xlabel('Teste');
ylabel('FOB');
yyaxis right
plot(aux, Resultado.final(aux,8));
ylabel('Tempo de simulação, s');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')

%% variaçãõ do eltismo
figure
aux = 14:16;
yyaxis left
plot(8:10, Resultado.final(aux,4), 'b' ), hold on, grid on
plot(8:10, Resultado.final(aux,7), '--r' )
title('Testes 8 a 10: Variação elitismo, [3, 10]');
xlabel('Teste');
ylabel('FOB');
yyaxis right
plot(8:10, Resultado.final(aux,8));
ylabel('Tempo de simulação, s');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')

SOLUCOES0 = unifrnd(3,10, 150, 3);
% Testes
N   = 40;        % N° de individuos
ger = 50;            % N° de iterações
MSG = 30;             % MaxStallGenerations
FT  = 1e-9;           % FunctionTolerance
EL  = ceil(0.1*N);   % EliteCount
CO  = 0.70;            % CrossoverFraction
Linf = 3;
Lsup =10;
% Estrutura de opções
SOLUCOES = SOLUCOES0(1:N, :);
opt = optimoptions(@ga, ...
    'PlotFcn', {@gaplotbestf},... @gaplotbestindiv}, ...
    'InitialPopulationMatrix', SOLUCOES,...
    'PopulationSize', N, ... (%%%)
    'MaxGenerations', ger, ... (%%%)
    'MaxStallGenerations', MSG, ... (%%%)
    'FunctionTolerance', FT, ... (%%%)
    'EliteCount', EL, ... (%%%)
    'CrossoverFraction', CO,... (%%%)
    'SelectionFcn', @selectionroulette,...
    'CrossoverFcn', @crossoversinglepoint, ...
    'MutationFcn', {@mutationadaptfeasible, 1/100}); % mutationuniform mutationadaptfeasible
% optimoptions gaoptimset     'Display', 'off',...
[x,fval, flag] = ga(@FOB_ga,d,[],[],[],[],Linf*[1 1 1],Lsup*[1 1 1],[],opt);

%% variaçãõ do crossver
figure
aux = 17:19;
aux2 = 11:13;
yyaxis left
plot(aux2, Resultado.final(aux,4), 'b' ), hold on, grid on
plot(aux2, Resultado.final(aux,7), '--r' )
title('Testes 11 a 13: Variação crossover, [3, 10]');
xlabel('Teste');
ylabel('FOB');
yyaxis right
plot(aux2, Resultado.final(aux,8));
ylabel('Tempo de simulação, s');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')

SOLUCOES0 = unifrnd(3,10, 150, 3);
% Testes
N   = 40;        % N° de individuos
ger = 50;            % N° de iterações
MSG = 30;             % MaxStallGenerations
FT  = 1e-9;           % FunctionTolerance
EL  = ceil(0.05*N);   % EliteCount
CO  = 0.90;            % CrossoverFraction
Linf = 3;
Lsup =10;
% Estrutura de opções
SOLUCOES = SOLUCOES0(1:N, :);
opt = optimoptions(@ga, ...
    'PlotFcn', {@gaplotbestf},... @gaplotbestindiv}, ...
    'InitialPopulationMatrix', SOLUCOES,...
    'PopulationSize', N, ... (%%%)
    'MaxGenerations', ger, ... (%%%)
    'MaxStallGenerations', MSG, ... (%%%)
    'FunctionTolerance', FT, ... (%%%)
    'EliteCount', EL, ... (%%%)
    'CrossoverFraction', CO,... (%%%)
    'SelectionFcn', @selectionroulette,...
    'CrossoverFcn', @crossoversinglepoint, ...
    'MutationFcn', {@mutationadaptfeasible, 1/100}); % mutationuniform mutationadaptfeasible
% optimoptions gaoptimset     'Display', 'off',...
[x,fval, flag] = ga(@FOB_ga,d,[],[],[],[],Linf*[1 1 1],Lsup*[1 1 1],[],opt);

%% CURVA
x1 = [-1	0	1	2	4	5	5	6];
x2 = [-2	-1	0	1	1	2	3	4];
y1 = [13	11	9	4	11	9	1	-1];

coef = [3.868428270879582   3.348360244117651  -6.242462566198178];
coef = [3.8828	3.3387	-6.2287];
coef = [3 3 3]

y2 = [ones(length(x1),1) x1' x2']*coef';

figure
plot3(x1, x2, y1, 'b'), hold on, grid on
plot3(x1, x2, y2, 'r')
legend('Curva original', 'Ajuste AG')
plot3(x1, x2, y1, 'b*')
plot3(x1, x2, y2, 'r*')

title('Curva do problema');
xlabel('A');
ylabel('B');
zlabel('C');














































