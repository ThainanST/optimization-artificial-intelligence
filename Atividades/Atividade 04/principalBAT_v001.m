% 02/05/2017
% T�cnicas Inteligentes - Bat Algorithm vc Algoritmo Gen�tico
% Prof. Ivo Chaves

clc
clear all
close all

%% Dados simula��o BAT
d = 3;          % Dimens�o
N = 40;         % N� de morcegos
iter = 50;     % N� de itera��es
alpha = 10^((log10(0.1))/(0.5*iter));  % (50%.iter , 0.1)
lambda = -log(1 - 0.9)/(0.3*iter);      % (30%.iter , 0.9)

fig = 0;
if fig == 1
    figure
    plot(1:iter, 1-exp(-lambda*(1:iter)) ), hold on, grid on
    plot(1:iter, cumprod(alpha*ones(1,iter)) )
    title(['r_i e A_i ( \lambda = ' num2str(lambda) ' , \alpha = ' num2str(alpha) ' )']);
    legend('Taxa de emiss�o, r_i', 'Amplitude sonora, A_i');
    xlabel('N�mero de itera��es');
    ylabel('Taxa de emiss�o e Amplitude');
end

% Inicializa��es
Linf = -10;
Lsup = 10;
SOLUCOES0 = unifrnd(Linf,Lsup, N, d);

%% Simula��es
% tic
% [x,fval] = batfunc( d, N, iter, y, a, SOLUCOES0, Linf, Lsup)
% toc
Nsimu = 50;

h = waitbar(0,'Buscando Solu��es... (0)');
for k = 1:Nsimu
    tic
    [x,fval] = batfunc( d, N, iter, lambda, alpha, SOLUCOES0, Linf, Lsup);
    BAT(k,1:d) = x;
    BAT(k,d+1) = fval;
    BAT(k,d+2) = toc;
    waitbar(k/Nsimu,h);
    ntext = ['Buscando Solu��o...' num2str(100*k/Nsimu, '(%0.f)')];
    set( get(findobj(h,'type','axes'),'title'), 'string', ntext);
    clear ntext
    clc
end

close(h)
[a, b] = sort(BAT(:,d+1));
BAT_sort_fval = BAT(b,:);

[a, b] = sort(BAT(:,end));
BAT_sort_time = BAT(b,:);

Resultado = [BAT_sort_fval([1 end],:);
            BAT_sort_time([1 end],:)];
Val_medio = mean([BAT(:,4) BAT(:,5)]);        
Resultado = [Resultado ones(4,1)*Val_medio(1) ones(4,1)*Val_medio(2)]
