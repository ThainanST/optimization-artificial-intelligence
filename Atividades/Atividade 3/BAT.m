% 23/04/2017
% Técnicas Inteligentes - Bat Algorithm
% Prof. Ivo Chaves
% Bat Algorithm

% Minimizar o erro de aproximação de uma curva
% Determinar por BAT os coeficientes da equação
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
% Possíveis valores de A, B e C
Dados_prob.A = [-10 10];
Dados_prob.B = Dados_prob.A;
Dados_prob.C = Dados_prob.A;
Dados_prob.d = 3;

% Tabela de dados de x1, x2 e yi
Dados_prob.x1x2y = [1  1 1 1 1 1 1 1; % A
                   -1  0 1 2 4 5 5 6; % B
                   -2 -1 0 1 1 2 3 4; % C
                    13 11 9 4 11 9 1 -1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Passo 1: Definição dos parâmetros do algoritmo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dados_simu.N = 15;       % N° de morcegos
Dados_simu.iter = 100;  % Critério de parada
tx = 0;
Dados_simu.a = 0.95;    % Tx amplitude sonora (0:1)
Dados_simu.y = 0.1;     % Tx emissão do pulso (0.001:0.9)


% Calculo de y em função do n° de iterações (alterar 0.9)
Dados_simu.y = -log(1 - 0.9)/(0.3*Dados_simu.iter);

% Calculo de alpha em função do n° de iterações (alterar 0.1)
Dados_simu.a = 10^((log10(0.1))/(0.5*Dados_simu.iter-2));

% Plot: Busca Local
% screenposition = get(fig,'Position')
screenposition = [59 145 1238 434];
fig = figure('Position', screenposition);

subplot(221)
plot(1:Dados_simu.iter, 1-exp(-Dados_simu.y*(1:Dados_simu.iter)) );
hold on
grid on
% Plot: Busca global
A(1) = 1;
for k = 2:Dados_simu.iter 
    A(k) = A(k-1)*Dados_simu.a;
end, clear k
plot(1:Dados_simu.iter, A ); hold on
title('r_i e A_i');
legend('Taxa de emissão, r_i', 'Amplitude sonora, A_i');

xlabel('Número de iterações');
ylabel('Taxa de emissão e Amplitude');
clear A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Passo 2: Inicialização dos morcegos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r = (b-a)*rand(1) + a
for k = 1:Dados_simu.N
    Solucao.MORCEGOS(k,:) = ...
    [(max(Dados_prob.A)-min(Dados_prob.A))*rand(1) + min(Dados_prob.A);...
    (max(Dados_prob.B)-min(Dados_prob.B))*rand(1) + min(Dados_prob.B);...
    (max(Dados_prob.C)-min(Dados_prob.C))*rand(1) + min(Dados_prob.C)];
end, clear k
% aux = Solucao.MORCEGOS;
% save('Morcegos0_30.mat', 'aux');
% clear aux

% load('Morcegos0_30')
% Solucao.MORCEGOS = aux; clear aux
% Solucao.MORCEGOS(1,:) = xbest;

% Limitações
for m = 1:Dados_simu.N
    Solucao.MORCEGOS(m,:) = limites(Solucao.MORCEGOS(m,:), min(Dados_prob.A), max(Dados_prob.A));
end, clear m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Passo 3: Cálculo da FOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m = 1:Dados_simu.N
    Solucao.MORCEGOS(m,Dados_prob.d+1) = FOB( Dados_prob, Solucao.MORCEGOS(m,1:Dados_prob.d) );
end, clear m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Passo 4: Atualização do melhor morcego
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[index, lin] = find(Solucao.MORCEGOS(:,end) == min(Solucao.MORCEGOS(:,end)) );
Solucao.MELHOR = Solucao.MORCEGOS(index, :);
clear lin index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Passo 5: Processo de iteração (parada pr iteração)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 1; % Iteração ZERO
% Inicialização: velocidade, tx emissao e ampliude para iteração 1
Solucao.r = zeros(1,Dados_simu.N);  % Tx de emissão (m,t,val)
Solucao.A = ones(1,Dados_simu.N);   % Amplitude sonora
Solucao.v = zeros(Dados_simu.N, Dados_prob.d);   % Velocidae Inicial

%% WAIT BAR
h = waitbar(0,'Buscando Solução... (0)');

while t <= Dados_simu.iter+tx
    % Laço secundário, morcegos
    for m = 1:Dados_simu.N
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Passo 7 a 9: Velocidade e deslocamento iniciais
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % atualização velocidade
    Solucao.v(m,:) = Solucao.v(m,:) + ...
        (Solucao.MELHOR(1:Dados_prob.d) - Solucao.MORCEGOS(m,1:Dados_prob.d) )*unifrnd(0,1);
    % Atualização deslocamento
    Solucao.TEMP(1:Dados_prob.d) = Solucao.MORCEGOS(m,1:Dados_prob.d) + Solucao.v(m,:);
    % Aplicação dos limites
    Solucao.TEMP(1:Dados_prob.d) = limites(Solucao.TEMP(1:Dados_prob.d), min(Dados_prob.A), max(Dados_prob.A));
    %% Passo 10 a 12:  Busca local
    if rand < Solucao.r(m)
        % Solucao.TEMP(1:Dados_prob.d) = Solucao.TEMP(1:Dados_prob.d) + 0.001*randn(1,Dados_prob.d);
        ep = unifrnd(-1, 1, 1, Dados_prob.d);
        Solucao.TEMP(1:Dados_prob.d) = Solucao.MELHOR(1:Dados_prob.d) + ep*mean(Solucao.A);
        % Aplicação dos limites
        Solucao.TEMP(1:Dados_prob.d) = limites(Solucao.TEMP(1:Dados_prob.d), min(Dados_prob.A), max(Dados_prob.A));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Passo 13: Perturbação em uma dimensão / Avaliação FOB
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rand < 0.8
    % sorteio dim
    Solucao.TEMP(1,randsample(1:Dados_prob.d,1)) =...
        unifrnd(min(Dados_prob.A), max(Dados_prob.A));
    end
    % Aplicação dos limites
    Solucao.TEMP = limites(Solucao.TEMP(1:Dados_prob.d), min(Dados_prob.A), max(Dados_prob.A));
    % Avaliação FOB
    Solucao.TEMP(1,Dados_prob.d +1) = FOB( Dados_prob, Solucao.TEMP );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Passo 14 a 18: Busca Global
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (rand < Solucao.A(m))||...
            (Solucao.TEMP(1,end)<=Solucao.MORCEGOS(m,end))
        Solucao.MORCEGOS(m,:) = Solucao.TEMP;       % Atualização do morcego
        Solucao.r(m) = 1-exp(-Dados_simu.y*t);      % atualização da taxa de emissao
        Solucao.A(m) = Dados_simu.a*Solucao.A(m);   % atualização da amplitude
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Passo 19: Atualiza MELHOR morcego
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    best = find(Solucao.MORCEGOS(:,end) == min(Solucao.MORCEGOS(:,end)));
    Solucao.MELHOR = Solucao.MORCEGOS(best(1),:);
    clear best
    end, clear m % end for
    
    % Armazenamento do melhor
    best = find(Solucao.MORCEGOS(:,end) == min(Solucao.MORCEGOS(:,end)));
    Solucao.BEST(t,:) = Solucao.MORCEGOS(best(1),:);
    clear best
    % Armazenamento do médio
    Solucao.MEAN(t,:) = mean(Solucao.MORCEGOS,1);
    % Armazenamentos todos morcegos
    Solucao.TODOS{t,1} = Solucao.MORCEGOS;
    
    % Plot3 dos dados
    subplot(2,2, [2 4])
    cla reset
    plot3(Solucao.TODOS{t,1}(:,1), Solucao.TODOS{t,1}(:,2), Solucao.TODOS{t,1}(:,3), 'k*'); hold on
    plot3(Solucao.BEST(t,1), Solucao.BEST(t,2), Solucao.BEST(t,3), 'r*')
    text(10,10,5, ...
        [ num2str(Dados_simu.N) ' Morcegos (' num2str(Solucao.BEST(t,1),'%.4f') ';'...
        num2str(Solucao.BEST(t,2),'%.4f') ';'...
        num2str(Solucao.BEST(t,3),'%.4f') ') - (' ...
        num2str(Solucao.BEST(t,4),'%.4f') ')'])
    grid on
    xlim(Dados_prob.A.*[0.8 1.2]);
    ylim(Dados_prob.A.*[0.8 1.2]);
    zlim(Dados_prob.A.*[0.8 1.2]);
    title(['Trajetória da solução, iteração ' num2str(t) ' / ' num2str(Dados_simu.iter)]);
    xlabel('Coef. A');
    ylabel('Coef. B');
    zlabel('Coef. C');
    
    subplot(221)
    plot(t, mean(Solucao.r), 'b*'); hold on, grid on
    plot(t, mean(Solucao.A), 'r*');
    title('Taxa de emissão e Amplitude');
    % legend('Taxa de emissão, r_i', 'Amplitude sonora, A_i');
    xlabel('Número de iterações');
    ylabel('r_i e A_i');
    
    xlim([0 Dados_simu.iter]);
    ylim([0 1])
    
    subplot(223)
    plot(t, Solucao.BEST(t,end),'b*'), hold on
    xlim([1 Dados_simu.iter]);
    title('Evolução da FOB')
    xlabel('Número de iterações');
    xlabel('FOB');
    
    drawnow
    F(t) = getframe(gcf);
    grid on
    
    %% WAIT BAR
    waitbar(t/Dados_simu.iter,h);
    ntext = ['Buscando Solução...' num2str(100*(t)/Dados_simu.iter, '(%0.f)')];
    set( get(findobj(h,'type','axes'),'title'), 'string', ntext);
    clear ntext
    t = t + 1;
end

video = VideoWriter('Trajetoria.avi', 'Uncompressed AVI');
video.FrameRate = 5;
open(video);
writeVideo(video,F);
close(video)

close(h);
clear h

Solucao.index = find(Solucao.MORCEGOS(:,end) == min(Solucao.MORCEGOS(:,end)) );

% close all
fprintf('*********************\n');
fprintf('Solução do problema\n');
fprintf('*********************\n');
fprintf('A = %.4f\nB = %.4f\nC = %.4f\n',Solucao.MORCEGOS(Solucao.index(1),1:3))
fprintf('*********************\n');
fprintf('Menor erro médio = %.3f\n',Solucao.MORCEGOS(Solucao.index(1),end));
fprintf('*********************\n');
fprintf('%d/%d morcego(s) chegou/chegaram à solução! \n',...
    length(Solucao.index), Dados_simu.N );
% Solucao.MORCEGOS

fprintf('Fim!!\n');


figure
plot3(Solucao.BEST(:,1), Solucao.BEST(:,2), Solucao.BEST(:,3), 'LineWidth', 2); hold on, grid on
% plot3(Solucao.MEAN(:,1), Solucao.MEAN(:,2), Solucao.MEAN(:,3), 'r');
legend('Melhor', 'Médio', 'Location','northwest');

plot3(Solucao.TODOS{1,1}(:,1), Solucao.TODOS{1,1}(:,2), Solucao.TODOS{1,1}(:,3), 'ko')
plot3(Solucao.TODOS{end,1}(:,1), Solucao.TODOS{end,1}(:,2), Solucao.TODOS{end,1}(:,3), 'k*')
plot3(Solucao.TODOS{end,1}(:,1), Solucao.TODOS{end,1}(:,2), Solucao.TODOS{end,1}(:,3), 'ko')

plot3([Solucao.TODOS{1,1}(:,1) Solucao.TODOS{100,1}(:,1)]', ...
      [Solucao.TODOS{1,1}(:,2) Solucao.TODOS{100,1}(:,2)]', ...
      [Solucao.TODOS{1,1}(:,3) Solucao.TODOS{100,1}(:,3)]', 'b--')

plot3(Solucao.BEST(1,1), Solucao.BEST(1,2), Solucao.BEST(1,3),'bo');
plot3(Solucao.BEST(end,1), Solucao.BEST(end,2), Solucao.BEST(end,3),'b*');

ntext = ['(' num2str(Solucao.MORCEGOS(Solucao.index(1),1))...
         ';' num2str(Solucao.MORCEGOS(Solucao.index(1),2))...
         ';' num2str(Solucao.MORCEGOS(Solucao.index(1),3)) ')'];
% text(Solucao.BEST(end,1), Solucao.BEST(end,2), Solucao.BEST(end,3), ntext); 
% plot3(Solucao.MEAN(1,1), Solucao.MEAN(1,2), Solucao.MEAN(1,3), 'ro')
% plot3(Solucao.MEAN(end,1), Solucao.MEAN(end,2), Solucao.MEAN(end,3), 'r*')

xlim(Dados_prob.A);
ylim(Dados_prob.A);
zlim(Dados_prob.A);

title(['Trajetória da solução' ntext]);
xlabel('Coef. A');
ylabel('Coef. B');
zlabel('Coef. C');
clear ntext




















