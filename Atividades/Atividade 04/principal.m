% 25/04/2017
% Técnicas Inteligentes - Bat Algorithm vc Algoritmo Genético
% Prof. Ivo Chaves
% Toolbox
% https://www.mathworks.com/help/gads/genetic-algorithm-options.html
% optimtool

clc
clear all
close all

%% DADOS BAT
d = 3;          % Dimensão
N = 15;         % N° de morcegos
iter = 200;     % N° de iterações
y = 0.0768;     % Lambda
a = 0.9532;     % Alpha

% Inicializações
Linf = -10;
Lsup = 10;
SOLUCOES0 = unifrnd(Linf,Lsup, N, d);
Nsimu = 10;

%% DADOS AG
opt = optimoptions(@ga,... 'PlotFcn', {@gaplotbestf @gaplotbestindiv}, ...     'FitnessScalingFcn', @fitscalingrank, ...
    'FitnessLimit', 1e-19,...
    'PopulationSize', N, ...
    'InitialPopulationMatrix', SOLUCOES0,...
    'SelectionFcn', @selectionroulette,...
    'EliteCount', ceil(0.3*N), ...
    'MaxGenerations', iter, ...
    'ConstraintTolerance', 1e-19);

h = waitbar(0,'Buscando Soluções... (0)');

%% SIMULAÇÕES
for k = 1:Nsimu
    tic
    [x,fval] = batfunc( d, N, iter, y, a, SOLUCOES0, Linf, Lsup);
    BAT(k,1:d) = x;
    BAT(k,d+1) = fval;
    BAT(k,d+2) = toc;
    
    tic
    [x,fval] = ga(@FOB_ga,d,[],[],[],[],Linf*[1 1 1],Lsup*[1 1 1], [], opt);
    clc
    AG(k,1:d) = x;
    AG(k,d+1) = fval;
    AG(k,d+2) = toc;
    
    waitbar(k/Nsimu,h);
    ntext = ['Buscando Solução...' num2str(100*k/Nsimu, '(%0.f)')];
    set( get(findobj(h,'type','axes'),'title'), 'string', ntext);
    clear ntext
end, clear k x fval

close(h)

[a, b] = sort(BAT(:,d+1));
BAT_sort = BAT(b,:);

[a, b] = sort(AG(:,d+1));
AG_sort = AG(b,:);

Linha = {'BAT - Maior';'BAT - Menor';'AG - Maior';'AG - Menor'};
TAB = [BAT(end,:); BAT(1,:); AG(end,:); AG(1,:)];

coefA = TAB(:,1);
coefB = TAB(:,2);
coefC = TAB(:,3);
coefFOB = TAB(:,4);
Tempo = TAB(:,5);

T = table(coefA,coefB,coefC,coefFOB,Tempo,'RowNames',Linha)
SOLUCOES0









