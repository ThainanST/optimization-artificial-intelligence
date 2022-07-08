%% ATIVIDADE_01_v005
% Ferramenta nnstart

clear all, close all, clc

var_simu.teste = 'IT31-0050';
var_simu.criar_rede = 1;
Treino.n_inf = 1;
var_simu.ft = 'trainlm'; 
% var_simu.ft = 'trainrp';
% var_simu.ft = 'trainlm';
% var_simu.ft = 'trainscg';
% var_simu.ft = 'trainbr';
var_simu.da = str2num(var_simu.teste([3])); % semanas anteriores
var_simu.dp = str2num(var_simu.teste([4])); % semanas posteriores
var_simu.nnco = str2num(var_simu.teste([6 7 8 9])); % Número de neurônios da camada oculta
fprintf('Teste: %s \nCriar rede = %d \nTreino = %s\n\n',...
    var_simu.teste,  var_simu.criar_rede, var_simu.ft);
%% ETAPA 01 - Leitura dos dados
dados.data = xlsread('\Dados\itapipoca2');
fprintf('ETAPA 01 (OK!) - Leitura dos dados \n')
% Database [70 30]
dados.horas_totais = floor(length(dados.data));
dados.dias_totais = floor(dados.horas_totais/24);
dados.semanas_totais = floor(dados.dias_totais/7);
dados.data_vals = [dados.horas_totais dados.dias_totais dados.semanas_totais];
dados.data0 = [dados.data(1,1) dados.data(1,2) dados.data(1,3)];
dados.dataend = [dados.data(end,1) dados.data(end,2) dados.data(end,3)];
clear dados.horas_totais dados.dias_totais dados.semanas_totais

dados.horas_train = floor(length(dados.data)*0.7);
dados.dias_train = floor(dados.horas_train/24);
dados.semanas_train = floor(dados.dias_train/7);
dados.train_vals = [dados.horas_train dados.dias_train dados.semanas_train];
dados.data_ini_train = [dados.data(1,1) dados.data(1,2) dados.data(1,3)];
dados.data_fim_train = [dados.data(dados.horas_train,1) dados.data(dados.horas_train,2) dados.data(dados.horas_train,3)];
clear dados.horas_train dados.dias_train dados.semanas_train

dados.horas_test = floor(length(dados.data)*0.3);
dados.dias_test = floor(dados.horas_test/24);
dados.semanas_test = floor(dados.dias_test/7);
dados.test_vals = [dados.horas_test dados.dias_test dados.semanas_test];
dados.data_ini_test = [dados.data(dados.train_vals(1)+24,1) dados.data(dados.train_vals(1)+24,2) dados.data(dados.train_vals(1)+24,3)];
dados.data_fim_test = [dados.data(end,1) dados.data(end,2) dados.data(end,3)];
clear dados.horas_test dados.dias_test dados.semanas_test

fprintf('BANCO DE DADOS (%d/%d/%d)-(%d/%d/%d)\n', dados.data0,dados.dataend);
fprintf('\nHoras = %d\nDias = %d\nSemanas = %d\n\n',...
    dados.data_vals(1), dados.data_vals(2), dados.data_vals(3));
% Treino (70%)
fprintf('TREINAMENTO(%d/%d/%d)-(%d/%d/%d)\n', dados.data_ini_train,dados.data_fim_train);
fprintf('\nHoras = %d\nDias = %d\nSemanas = %d\n\n',...
    dados.train_vals(1), dados.train_vals(2), dados.train_vals(3));
% Teste (30%)
fprintf('TESTE (%d/%d/%d)-(%d/%d/%d)\n', dados.data_ini_test,dados.data_fim_test);
fprintf('\nHoras = %d\nDias = %d\nSemanas = %d\n\n',...
    dados.test_vals(1), dados.test_vals(2), dados.test_vals(3));
%% ETAPA 02 - Dados de treinamento e Estratégia de rede
% 16/04/2016 - 07/04/2017
% Y	 M	D	H	T	P	W
Treino.posicao0 = find( dados.data(1:end,1) == dados.data0(1) & ... % DIA
                 dados.data(1:end,2) == dados.data0(2) & ... % MES
                 dados.data(1:end,3) == dados.data0(3) & ... % ANO
                 dados.data(1:end,4) == 0);            % HORA

% N° de treinamentos
Treino.n_trains = dados.train_vals(3) - var_simu.da + 1;
% N° de horas de um treinamento
Treino.horas_train = var_simu.da*24;
% N° de intervalo de Semanas
Treino.n_inter = 3*24; % 3 dias, 24 horas

% vetor de posições iniciais de cada hora [TEMP PRESS WIND]
Treino.pos_aux = 1:Treino.n_inf:Treino.n_inf*Treino.horas_train;
% posição inicial e final do treinamento
Treino.ini = Treino.posicao0;
Treino.fim = Treino.horas_train;

Treino = FORM_in_out( Treino, var_simu, dados);

close all
dias = 8;
figure
plot(dados.data(1:dias*24,7))
title(['Perfil de Vento de ' num2str(dias) ' dias'])
xlabel('Horas'), ylabel('Vento')

figure
cor{1} = '-bs';
cor{2} = '-rs';
cor{3} = '-ms';

z = zeros(size(Treino.S(:,1),1),1);
z1 = [z; z; z];
plot(Treino.E(:,1), cor{1}), hold on
plot([z1 ; Treino.S(:,1)], cor{1},'LineWidth', 2)

z2 = [z; z; z; z];
plot([z2; Treino.E(:,2)], cor{2}), hold on
plot([[z2; z1] ; Treino.S(:,2)], cor{2},'LineWidth', 2)
title('Treino 1')
xlabel('Horas'), ylabel('Vento')

%% Criação de dados de teste
% posiçã inicial dos testes
Teste.posicao0 = dados.train_vals(1);
% N° de testes
Teste.n_trains = 15;
% N° de horas de um teste
Teste.horas_teste = var_simu.da*24;
% N° de intervalo
Teste.n_inter = 3*24; % 3 dias, 24 horas
% vetor de posições iniciais de cada hora [TEMP PRESS WIND]
Teste.n_inf = Treino.n_inf;
Teste.pos_aux = 1 : Teste.n_inf : Teste.n_inf*Teste.horas_teste;
% posição inicial e final do treinamento
Teste.ini = Teste.posicao0;
Teste.fim = Teste.posicao0 + Teste.horas_teste - 1;

% Laço de criação da entrada
Teste = FORM_in_out( Teste, var_simu, dados);

save('Teste.mat', 'Teste');
E0 = Treino.E;
S0 = Treino.S;

E = Teste.E(:,1:10);
S = Teste.S(:,1:10);

%% ETAPA 03 - Criação e treinamento da rede

if var_simu.criar_rede == 1
    
    net = feedforwardnet(var_simu.nnco, var_simu.ft); % Definição da rede
    net.layers{1}.transferFcn = 'purelin';  % Função de ativação da camada de entrada
    net.layers{2}.transferFcn = 'purelin';  % Função de ativação da camada de saída
    net.trainParam.epochs   = 2*1e3;     % NUMERO DE EPOCAS
    net.trainParam.goal     = 1*1e-6;    % ERRO?
    net.trainParam.min_grad = 1*1e-9;    % GRADIENTE DO ERRO
    net = initlay(net); % net = init(net);
    
    % net.trainParam.showWindow = false; % DESABILITAR NNTRAINTOOL
    [net,tr] = train(net,Treino.E,Treino.S);
    save([pwd '\Resultados - INMET\' var_simu.teste '.mat'],'net', 'tr');
else
    load([pwd '\Resultados - INMET\' var_simu.teste '.mat']);
end

fprintf('ETAPA 03 (OK!) - Criação e treinamento da rede \n')
%% ETAPA 04 - Teste da Rede

% Teste
for k = 1:4 %size(Testes.E,1)
    Teste.s(:,k) = sim(net,Teste.E(k,1));
end

% save('Testes.mat', 'Testes')
% E1 = Testes.E;
% S1 = Testes.S;






n_plot = 4;

for k = 1:4
    subplot(n_plot/2,2,k);
    v_real = Teste.S(:,k); % coluna vento (17)
    v_prev = Teste.s(:,k);                % Vento previsto
    
    % Aplicação dos critérios (NIAE correlação erro)
    %[erro_p(:,k), NIAE_val(k,:) ] = NIAE_function([ha 23], h(1,ha+1:end), v_real(1,ha+1:end), v_prev);
    
    % Plotagem
    subplot(211)
    plot(v_real, 'b-o', 'Linewidth', 2); grid on; hold on
    subplot(212)
    plot(v_prev, 'r-o', 'Linewidth', 2);
    xlim([0 23]), ylim([0 10]);
    legend('Real','Previsto');
    xlabel('Tempo (hora)');
    ylabel('Velocidade (m/s)')    
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    






















% % Dias do teste (sempre par) 
% % dia do teste(2 04 2016) + 8meses
% datas = [ 10 01 2017
%           20 01 2017
%           10 02 2017
%           20 02 2017];
% 
% [n_plot, aux] = size(datas); % Número de plots (sempre par)
% posicao = zeros(n_plot, 24); % Inicialização
% 
% % Posição dos dados de vento nos dados
% for k = 1:n_plot
%     posicao(k,:) = find( dados(1:end,1) == datas(k,3) & ... % ANO
%                          dados(1:end,2) == datas(k,2) & ... % MES
%                          dados(1:end,3) == datas(k,1) )';   % DIA
% end, clear k
% 
% % Configuração das entradas e aplicação da rede
% for m = 1:n_plot
%     % Configuração das entradas
%     for k = 1:length(posicao)-ha
%         ini = posicao(m,1)-1+k;
%         fim = ini + ha - 1;
%         E(:,k) = dados(ini:fim,17); % coluna vento (17)
%     end
%     % Aplicação da rede
%     aux = sim(net,E);
%     for n = 1:size(aux, 1)
%         S(m,:,n) = aux(n,:);
%     end
%     clear k n ini fim aux
% end
% 
% fprintf('ETAPA 04 (OK!) - Teste da Rede \n')













%%
% f = 60;
% t = 0:1/(10*f):1/f;
% w = 2*pi*60;
% y = sin(w*t);
% y1 = sin(w*t+pi/6);
% err = y-y1;
% 
% figure(1)
% subplot(211), plot(t, y), hold on, plot(t, y1, 'r'), xlim([t(1) t(end)]),ylim([-2 2])
% subplot(212)
% E = errorbar(t,y,err,'-o', 'MarkerSize',4,...
%     'MarkerEdgeColor','blue','MarkerFaceColor','blue')
% xlim([t(1) t(end)]), ylim([-2 2])
