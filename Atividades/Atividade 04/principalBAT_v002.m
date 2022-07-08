% 02/05/2017
% Técnicas Inteligentes - Bat Algorithm vc Algoritmo Genético
% Prof. Ivo Chaves

clc
clear all
close all
% iter	N	?	?
teste = [...
50	40	0.50	0.30
100	40	0.50	0.30
150	40	0.50	0.30
200	40	0.50	0.30
100	10	0.50	0.30
100	40	0.50	0.30
100	70	0.50	0.30
100	100	0.50	0.30
100	40	0.25	0.30
100	40	0.50	0.30
100	40	0.75	0.30
100	40	0.50	0.15
100	40	0.50	0.30
100	40	0.50	0.45];

Linf = 3;
Lsup = 10;
d = 3;             % Dimensão
Nsimu = 25;

% Inicializações
SOLUCOES0 = unifrnd(Linf,Lsup, 100, d);

h = waitbar(0,'Buscando Soluções... (0)');
for t = 1:size(teste,1)
    %% Dados simulação BAT
    N = teste(t,2);    % N° de morcegos
    iter = teste(t,1); % N° de iterações
    alpha = 10^((log10(0.1))/(teste(t,3)*iter));  % (50%.iter , 0.1)
    lambda = -log(1 - 0.9)/(teste(t,4)*iter);     % (30%.iter , 0.9)   

    
    %% Simulações
    solucoes0 = SOLUCOES0(1:N,:);
    for k = 1:Nsimu
        tic
        [x,fval, FOB] = batfunc( d, N, iter, lambda, alpha, solucoes0, Linf, Lsup);
        BAT(k,1:d) = x;
        BAT(k,d+1) = fval;
        BAT(k,d+2) = toc;
    end
    
    [a, b] = sort(BAT(:,d+1));
    BAT_sort_fval = BAT(b,:);
    clear a b
    
    [a, b] = sort(BAT(:,end));
    BAT_sort_time = BAT(b,:);
    clear a b
    
    aux = [BAT_sort_fval([1 end],:);
        BAT_sort_time([1 end],:)];
    Val_medio = mean([BAT(:,4) BAT(:,5)]);
    Resultado.val{t,1} = [aux ones(4,1)*Val_medio(1) ones(4,1)*Val_medio(2)];
    Resultado.teste = teste;
    
    clear x fval BAT BAT_sort_fval BAT_sort_time aux Val_medio
     
    waitbar(t/length(teste),h);
    ntext = ['Buscando Solução...' num2str(100*t/length(teste), '(%0.f)')];
    set( get(findobj(h,'type','axes'),'title'), 'string', ntext);
    clear ntext
    
end

Resultado.SOL0 = SOLUCOES0;
close(h)

for k = 1:length(teste)
Resultado.final(k,:) = Resultado.val{k}(1,:);
end


%% Testes 1,4 e 5,8
figure
subplot(121)
aux = 1:4;
yyaxis left
plot(aux, Resultado.final(aux,4), 'b' ), hold on, grid on
plot(aux, Resultado.final(aux,6), '--k' )
ylabel('FOB');
xlabel('Teste');
ylim([1.9218 1.9219])
yyaxis right
plot(aux, Resultado.final(aux,7));
ylabel('Tempo de simulação, s');
title('Testes 1 a 4: Aumento iterações, [3, 10]');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')


subplot(122)
aux = 5:8;
yyaxis left
plot(aux, Resultado.final(aux,4), 'b' ), hold on, grid on
plot(aux, Resultado.final(aux,6), '--r' )
title('Testes 5 a 8: Aumento Morcegos, [3, 10]');
xlabel('Teste');
ylabel('FOB');
yyaxis right
plot(aux, Resultado.final(aux,7));
ylabel('Tempo de simulação, s');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')

%% Evolução da FOB (iterações)
SOLUCOES0 = unifrnd(3,10, 40, 3);
iter = [50 100 150 200];
alpha = 10.^((log10(0.1))./(0.5*iter));
lambda = -log(1 - 0.9)./(0.3*iter);

figure
for m = 1:4
    subplot(2,2,m)
    for k=1:5
        [x,fval, FOB] = batfunc( 3, 40, iter(m), lambda(m), alpha(m), SOLUCOES0, 3, 10);
        plot(FOB(:,4)), hold on, grid on
    end
    title(['Evolução da FOB (' num2str(iter(m)) ' iterações)']);
    xlabel('Número de iterações');
    ylabel('FOB');
end

%% Evolução da FOB Morcegos
M = [10 40 70 100];
iter = 100;
alpha = 10.^((log10(0.1))./(0.5*iter));
lambda = -log(1 - 0.9)./(0.3*iter);

figure
for m = 1:4
    subplot(2,2,m)
    SOLUCOES0 = unifrnd(3,10, M(m), 3);
    for k=1:5
        [x,fval, FOB] = batfunc( 3, M(m), 100, 0.0767, 0.95499, SOLUCOES0, 3, 10);
        plot(FOB(:,4)), hold on, grid on
    end
    title(['Evolução da FOB (' num2str(M(m)) ' morcegos)']);
    xlabel('Número de iterações');
    ylabel('FOB');
end

%% Amplitude
figure
subplot(121)
alpha = 10.^((log10(0.1))./([0.25 0.5 0.75]*100));
% 0.9120    0.9550    0.9698
plot(1:100, cumprod(alpha(1)*ones(1,100)) ), hold on, grid on
plot(1:100, cumprod(alpha(2)*ones(1,100)),'r' )
plot(1:100, cumprod(alpha(3)*ones(1,100)),'m' )
title('Variação da taxa de drecréscimo de amplitude');
legend('Teste 9, \alpha (25%)', 'Teste 10, \alpha (50%)', 'Teste 11, \alpha (75%)');
xlabel('Número de iterações');
ylabel('Amplitude da onda sonora');

subplot(122)
aux = 9:11;
yyaxis left
plot(aux, Resultado.final(aux,4), 'b' ), hold on, grid on
plot(aux, Resultado.final(aux,6), '--r' )
title('Testes 9 a 11: Aumento de \alpha, [3, 10]');
xlabel('Teste');
ylabel('FOB');
yyaxis right
plot(aux, Resultado.final(aux,7));
ylabel('Tempo de simulação, s');
legend('FOB', 'FOB Média do teste','Tempo médio simulação')

figure
SOLUCOES0 = unifrnd(3,10, 40, 3);
figure
for m = 1:3
    subplot(2,2,m)
    for k=1:5
        [x,fval, FOB] = batfunc( 3, 40, 100, 0.0768, alpha(m), SOLUCOES0, 3, 10);
        plot(FOB(:,4)), hold on, grid on
    end
    title(['Evolução da FOB ( \alpha = ' num2str(alpha(m)) ')']);
    xlabel('Número de iterações');
    ylabel('FOB');
end

figure
for m = 1:3
    for k=1:5
        [x,fval, FOB, Am, rm] = batfunc( 3, 40, 100, 0.0768, alpha(m), SOLUCOES0, -10, 10);
        Am_0(:,k) = Am;
        rm_0(:,k) = rm;
    end
    plot(mean(Am_0,2)), hold on, grid on
%     plot(mean(rm_0,2)), hold on, grid on
    title(['Evolução da FOB ( \alpha = ' num2str(alpha(m)) ')']);
    xlabel('Iterações');
    ylabel('Amplitude da onda');
end



%% Emissão de pulso
subplot(121)
lambda = -log(1 - 0.9)./([0.15 0.30 0.45]*100);
% 0.1535    0.0768    0.0512
plot(1:100, 1-exp(-lambda(1)*(1:100)) ), hold on, grid on
plot(1:100, 1-exp(-lambda(2)*(1:100)), 'r' )
plot(1:100, 1-exp(-lambda(3)*(1:100)), 'm' )
title('Variação da taxa de emissão de pulso');
legend('Teste 12, \lambda (15%)', 'Teste 13, \lambda (30%)', 'Teste 14, \lambda (45%)');
xlabel('Número de iterações');
ylabel('Taxa de emissão de pulso');

subplot(122)
aux = 12:14;
yyaxis left
plot(aux, Resultado.final(aux,4), 'b' ), hold on, grid on
plot(aux, Resultado.final(aux,6), '--r' )
title('Testes 12 a 14: Aumento de \lambda, [3, 10]');
xlabel('Teste');
ylabel('FOB');
yyaxis right
plot(aux, Resultado.final(aux,7));
ylabel('Tempo de simulação, s');
legend('FOB', 'FOB Média do teste','Tempo simulação')

figure
SOLUCOES0 = unifrnd(3,10, 40, 3);
for m = 1:3
    subplot(2,2,m)
    for k=1:5
        [x,fval, FOB] = batfunc( 3, 40, 100, lambda(m), 0.9550, SOLUCOES0, 3, 10);
        plot(FOB(:,4)), hold on, grid on
    end
    title(['Evolução da FOB ( \lambda = ' num2str(lambda(m)) ')']);
    xlabel('Número de iterações');
    ylabel('FOB');
end

%% CURVA
x1 = [-1	0	1	2	4	5	5	6];
x2 = [-2	-1	0	1	1	2	3	4];
y1 = [13	11	9	4	11	9	1	-1];

coef = [3.868428270879582   3.348360244117651  -6.242462566198178];
coef = [3 3 3]

y2 = [ones(length(x1),1) x1' x2']*coef';

figure
plot3(x1, x2, y1, 'b'), hold on, grid on
plot3(x1, x2, y2, 'r')
legend('Curva original', 'Ajuste BAT')
plot3(x1, x2, y1, 'b*')
plot3(x1, x2, y2, 'r*')

title('Curva do problema');
xlabel('A');
ylabel('B');
zlabel('C');





























