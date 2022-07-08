% 02/05/2017
% Técnicas Inteligentes - Bat Algorithm vc Algoritmo Genético
% Prof. Ivo Chaves

clc
clear all
close all

time = 0;
Nsimu = 25;  % N° de simulações
d = 3;       % N° de dimensões
Linf = -10;  % Limite inferior 3;%
Lsup = 10;   % Limite superior
SOLUCOES0 = unifrnd(Linf,Lsup, 150, d);
% Testes
N   = 40;           % N° de individuos
ger = 50;           % N° de iterações
MSG = 30;           % MaxStallGenerations
FT  = 1e-9;         % FunctionTolerance
EL  = ceil(0.05*N); % EliteCount
CO  = 0.7;          % CrossoverFraction

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
    'CrossoverFcn', @crossovertwopoint, ...
    'MutationFcn', {@mutationuniform, 1/100} ); % mutationuniform mutationadaptfeasible
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
Resultado.val{1} = [Resul ones(4,1)*Val_medio(1) ones(4,1)*Val_medio(2)];
clc

for k = 1:1
    Resultado.final(k,:) = Resultado.val{k}(1,:)
end

Resultado.final

resultado = [
3.9243	3.3227	-6.2152	0.0748	0.1600	0	0.0778	0.1908
3.9865	3.3055	-6.2049	0.0750	0.1621	0	0.0987	0.1890
3.2147	3.9816	-7.0643	0.0807	0.1552	0	0.1688	0.1880
3.9002	3.3358	-6.2291	0.0748	0.1703	0	0.0996	0.2030
3.3222	3.7663	-6.7300	0.0770	0.3665	0	0.1010	0.2950
3.9582	3.3125	-6.2084	0.0749	0.1782	0	0.0990	0.1904
3.6268	3.5290	-6.4512	0.0757	0.1793	0	0.0844	0.2123
3.9375	3.3186	-6.2124	0.0748	0.1783	0	0.1055	0.2207
3.8788	3.3465	-6.2397	0.0748	0.2389	0	0.0772	0.2761
];



%% Testes 1,5
figure
aux = 1:5;
yyaxis left
plot(aux, resultado(aux,4), 'b' ), hold on, grid on
plot(aux, resultado(aux,7), '--k' )
ylabel('FOB');
xlabel('Teste');
yyaxis right
plot(aux, resultado(aux,8));
ylabel('Tempo de simulação, s');
title('Testes 1 a 5: Variação tipo seleção, [-10, 10]');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')

%% Testes 1,4
figure
aux = 6:9;
aux2 = 1:4;
yyaxis left
plot(aux2, resultado(aux,4), 'b' ), hold on, grid on
plot(aux2, resultado(aux,7), '--k' )
ylabel('FOB');
xlabel('Teste');
yyaxis right
plot(aux2, resultado(aux,8));
ylabel('Tempo de simulação, s');
title('Testes 1 a 4: Variação tipo cruzamento, [-10, 10]');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')














































