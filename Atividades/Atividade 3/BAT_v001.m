% 23/04/2017
% T�cnicas Inteligentes - Bat Algorithm
% Prof. Ivo Chaves
% Bat Algorithm

% Minimizar o erro de aproxima��o de uma curva
% Determinar por BAT os coeficientes da equa��o
% y = A + B.x1 + C.x2
% x1 = [-1  0 1 2 4  5 5  6]
% x2 = [-2 -1 0 1 1  2 3  4]
% yi = [13 11 9 4 11 9 1 -1];
% com: 3 <= A, B e C <= 10

%% Dados do problema
clear all
close all
clc
% xbest = [4.3927 3.3583 -6.4782];
% Poss�veis valores de A, B e C
A = [3 10];
B = A;
C = A;
d = 3;

% Tabela de dados de x1, x2 e yi
x1x2y = [-1 0 1 2 4 5 5 6;...
    -2 -1 0 1 1 2 3 4; 13 11 9 4 11 9 1 -1];

%% Passo 1: Defini��o dos par�metros do algoritmo
N = 30;       % N� de morcegos
iter = 100;  % Crit�rio de parada
a = 0.95;    % Tx amplitude sonora (0:1)
y = 0.1;     % Tx emiss�o do pulso (0.001:0.9)


% Calculo de y em fun��o do n� de itera��es (alterar 0.9)
y = -log(1 - 0.9)/(0.3*iter);

% Calculo de alpha em fun��o do n� de itera��es (alterar 0.1)
a = 10^((log10(0.25))/(0.5*iter-2));

% Plot: Busca Local
% screenposition = get(fig,'Position')
screenposition = [277 273 1238 434];
fig = figure('Position', screenposition);

subplot(121)
plot(1:iter, 1-exp(-y*(1:iter)) );
hold on
grid on
% Plot: Busca global
Av(1) = 1;
for k = 2:iter 
    Av(k) = Av(k-1)*a;
end, clear k
plot(1:iter, Av ); hold on
title('r_i e A_i');
legend('Taxa de emiss�o, r_i', 'Amplitude sonora, A_i');

xlabel('N�mero de itera��es');
ylabel('Taxa de emiss�o e Amplitude');
clear Av
%% Passo 2: Inicializa��o dos morcegos
% r = (b-a)*rand(1) + a
% for k = 1:N
%     MORCEGOS(k,:) = ...
%     [(max(A)-min(A))*rand(1) + min(A);...
%     (max(B)-min(B))*rand(1) + min(B);...
%     (max(C)-min(C))*rand(1) + min(C)];
% end, clear k
% aux = MORCEGOS;
% save('Morcegos0_30.mat', 'aux');
% clear aux

load('Morcegos0_30')
MORCEGOS = aux; clear aux
% MORCEGOS(1,:) = xbest;

for m = 1:N
    MORCEGOS(m,:) = limites(MORCEGOS(m,:), min(A), max(A));
end

%% Passo 3: Avalia��o do melhor morcego
MORCEGOS = FOBinicial_v001( N, d, x1x2y, MORCEGOS );

%% Passo 4: Atualiza��o do melhor morcego
[index, lin] = find(MORCEGOS(:,end) == ...
    min(MORCEGOS(:,end)));
MELHOR = MORCEGOS(index, :);
clear lin index

%% Passo 5: Processo de itera��o (parada pr itera��o)

t = 1; % Itera��o ZERO
% Inicializa��o: velocidade, tx emissao e ampliude para itera��o 1
r = zeros(1,N);  % Tx de emiss�o (m,t,val)
Amp = ones(1,N);   % Amplitude sonora
v = zeros(N, d);   % Velocidae Inicial

%% WAIT BAR
h = waitbar(0,'Buscando Solu��o... (0)');

while t <= iter
    % La�o secund�rio, morcegos
    for m = 1:N
    %% Passo 7 a 9: Velocidade e deslocamento iniciais
    % atualiza��o velocidade
    v(m,:) = v(m,:) + (MELHOR(1:d) - MORCEGOS(m,1:d) )*rand(1);
    % Atualiza��o deslocamento
    TEMP(1:d) = MORCEGOS(m,1:d) + v(m,:);
    % Aplica��o dos limites
    TEMP(1:d) = limites(TEMP(1:d), min(A), max(A));
    %% Passo 10 a 12:  Busca local
    if rand < r(m)
%         TEMP(1:d) = TEMP(1:d) + 0.001*randn(1,d);
        e = (1-(-1))*rand(1) + (-1);
        TEMP(1:d) = TEMP(1:d) + e*mean(Amp);
        clear e
        % Aplica��o dos limites
        TEMP(1:d) = limites(TEMP(1:d), min(A), max(A));
    end
    %% Passo 13: Perturba��o em uma dimens�o / Avalia��o FOB
    % sorteio dim
    dim = randsample(1:d, 1);
    TEMP(1,dim) = (0.5+0.5)*rand(1) -0.5; clear dim
    % Avalia��o FOB
    fob = FOB_v001( x1x2y, d, TEMP );
    TEMP(1,d +1) = fob;
    clear fob
    %% Passo 14 a 18: Busca Global
    if (Amp(m) > rand) | ( TEMP(1,d+1) <= MORCEGOS(m,d+1) )
        MORCEGOS(m,:) = TEMP(1,:);
        % atualiza��o da taxa de emissao
        r(m) = 1-exp(-y*t);
        % atualiza��o da amplitude
        Amp(m) = a*Amp(m);
    end
    %% Passo 19: Atualiza MELHOR
    best = find(MORCEGOS(:,end) == min(MORCEGOS(:,end)));
    MELHOR(t,:) = MORCEGOS(best,:);
    clear best
    end, clear m
    % Armazenamento do melhor da itera��o
    best = find(MORCEGOS(:,end) == min(MORCEGOS(:,end)));
    BEST(t,:) = MORCEGOS(best,:);
    clear best
    % Armazenamento do m�dio
    MEAN(t,:) = mean(MORCEGOS,1);
    % Armazenamentos todos morcegos
    TODOS{t,1} = MORCEGOS;
    
    % Plot3 dos dados
    subplot(122)
    plot3(TODOS{t,1}(:,1), TODOS{t,1}(:,2), TODOS{t,1}(:,3), 'k*')
    grid on
    xlim(A);
    ylim(A);
    zlim(A);
    title(['Trajet�ria da solu��o, itera��o ' num2str(t)]);
    xlabel('Coef. A');
    ylabel('Coef. B');
    zlabel('Coef. C');
   
    
    subplot(121)
    plot(t, mean(r), 'b*'); hold on, grid on
    plot(t, mean(Amp), 'r*');
    title('r_i e A_i');
    % legend('Taxa de emiss�o, r_i', 'Amplitude sonora, A_i');
    xlabel('N�mero de itera��es');
    ylabel('Taxa de emiss�o e Amplitude');
    
    xlim([0 iter]);
    ylim([0 1])
    
    drawnow
    F(t) = getframe(gcf);
    pause(0.01)
    grid on
    
    %% WAIT BAR
    waitbar(t/iter,h);
    ntext = ['Buscando Solu��o...' num2str(100*(t)/iter, '(%0.f)')];
    set( get(findobj(h,'type','axes'),'title'), 'string', ntext);
    clear ntext
    t = t + 1;
end

video = VideoWriter('Trajetoria.avi', 'Uncompressed AVI');
open(video);
writeVideo(video,F);
close(video)

close(h);
clear h

index = find(MORCEGOS(:,end) == min(MORCEGOS(:,end)) );

% close all
fprintf('*********************\n');
fprintf('Solu��o do problema\n');
fprintf('*********************\n');
fprintf('A = %d\nB = %d\nC = %d\n',MORCEGOS(index(1),1:3))
fprintf('*********************\n');
fprintf('Menor erro m�dio = %.3f\n',MORCEGOS(index(1),end));
fprintf('*********************\n');
fprintf('%d/%d morcego(s) chegou/chegaram � solu��o! \n',...
    length(index), N );
% MORCEGOS

fprintf('Fim!!\n');


figure
plot3(BEST(:,1), BEST(:,2), BEST(:,3)); hold on, grid on
plot3(MEAN(:,1), MEAN(:,2), MEAN(:,3), 'r');
legend('Melhor', 'M�dio', 'Location','northwest');

plot3(TODOS{1,1}(:,1), TODOS{1,1}(:,2), TODOS{1,1}(:,3), 'ko')
plot3(TODOS{end,1}(:,1), TODOS{end,1}(:,2), TODOS{end,1}(:,3), 'm*')

plot3(BEST(1,1), BEST(1,2), BEST(1,3),'bo');
plot3(BEST(end,1), BEST(end,2), BEST(end,3),'b*');
ntext = ['(' num2str(MORCEGOS(index(1),1))...
         ';' num2str(MORCEGOS(index(1),2))...
         ';' num2str(MORCEGOS(index(1),3)) ')'];
% text(BEST(end,1), BEST(end,2), BEST(end,3), ntext); 
plot3(MEAN(1,1), MEAN(1,2), MEAN(1,3), 'ro')
plot3(MEAN(end,1), MEAN(end,2), MEAN(end,3), 'r*')

xlim([2.5 7])
ylim([2.5 7])
zlim([2.5 7])

title(['Trajet�ria da solu��o' ntext]);
xlabel('Coef. A');
ylabel('Coef. B');
zlabel('Coef. C');
clear ntext




















