% 02/05/2017
% Técnicas Inteligentes - Bat Algorithm vc Algoritmo Genético
% Prof. Ivo Chaves

clc
clear all
close all

%% Testes
N = 75;      % N° de individuos
ger = 100;   % N° de iterações
time = 0;
MSG = 25;
FT = 1e-9;
EL = ceil(0.05*N);
CO = 0.9;

Nsimu = 25;  % N° de simulações
solInit = 1; % 1 = criar; 0 carregar

%% Dados
d = 3;       % N° de dimensões
Linf = -10;  % Limite inferior 3;%
Lsup = 10;   % Limite superior
% Inicializações
name = ['Morcegos0_' num2str(N) '.mat'];
if solInit == 1
    % Criar soluções
    SOLUCOES = unifrnd(Linf,Lsup, N, d);
    save(name, 'SOLUCOES');
else
    % Carregar soluções
    load(name);
end

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

h = waitbar(0,'Buscando Soluções... (0)');
for k = 1:Nsimu
    tic
    [x,fval, flag, output] = ga(@FOB_ga,d,[],[],[],[],Linf*[1 1 1],Lsup*[1 1 1],[],opt);
    AG(k,1:d) = x;
    AG(k,d+1) = fval;
    AG(k,d+2) = toc;
    AG(k,d+3) = flag;
    waitbar(k/Nsimu,h);
    ntext = ['Buscando Solução...' num2str(100*k/Nsimu, '(%0.f)')];
    set( get(findobj(h,'type','axes'),'title'), 'string', ntext);
    clear ntext
    clc
end

close(h)
[a, b] = sort(AG(:,d+1));
AG_sort_fval = AG(b,:);

[a, b] = sort(AG(:,end));
AG_sort_time = AG(b,:);

Resultado = [AG_sort_fval([1 end],:);
            AG_sort_time([1 end],:)];
Val_medio = mean([AG(:,4) AG(:,5)])        
Resultado = [Resultado ones(4,1)*Val_medio(1) ones(4,1)*Val_medio(2)]
        

    


















