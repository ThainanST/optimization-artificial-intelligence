%% ATIVIDADE_01_v005
% Ferramenta nnstart

clear all, close all, clc

teste = 'IT71-0100'; % A
criar_rede = 1;
ft = 'trainrp';
% ft = 'trainlm';
% ft = 'trainscg';
% ft = 'trainbr';
fprintf('Teste: %s \nCriar rede = %d \nTreino = %s\n\n', teste,  criar_rede, ft)
%% INFORMAÇÔES

% Tipo de treino
        % trainlm is a network training function that updates weight and
        %       bias states according to Levenberg-Marquardt optimization.
        % trainrp is a network training function that updates weight and bias
        %       values according to the resilient backpropagation algorithm (RPROP).
        % trainbr trains a network with weight and bias learning rules
        %       with batch updates. The weights and biases are updated at the end of
        %       an entire pass through the input data.  It
        %       minimizes a combination of squared errors and weights
        %       and, then determines the correct combination so as to produce a
        %       network which generalizes well.  The process is called Bayesian
        %       regularization
        % trainscg is a network training function that updates weight and
                % bias values according to the scaled conjugate gradient method.

    %      tansig Symmetric sigmoid transfer function
    %      purelin Linear transfer function
    % 'tansig' 'purelin'               
    % net.performFcn = 'mse';  %Name of a network performance function %type help nnperformance
    % view(net)
%     net = init(net)
%     net = initzero(net, [1 1]); % initzero(S,PR)     
%% ETAPA 01 - Leitura dos dados
dados = xlsread('\Dados\itapipoca2');
fprintf('ETAPA 01 (OK!) - Leitura dos dados \n')
% Database [70 30]
horas_totais = floor(length(dados));
dias_totais = floor(horas_totais/24);
semanas_totais = floor(dias_totais/7);
data_vals = [horas_totais dias_totais semanas_totais];
data0 = [dados(1,1) dados(1,2) dados(1,3)];
dataend = [dados(end,1) dados(end,2) dados(end,3)];
clear horas_totais dias_totais semanas_totais

horas_train = floor(length(dados)*0.7);
dias_train = floor(horas_train/24);
semanas_train = floor(dias_train/7);
train_vals = [horas_train dias_train semanas_train];
data_ini_train = [dados(1,1) dados(1,2) dados(1,3)];
data_fim_train = [dados(horas_train,1) dados(horas_train,2) dados(horas_train,3)];
clear horas_train dias_train semanas_train

horas_test = floor(length(dados)*0.3);
dias_test = floor(horas_test/24);
semanas_test = floor(dias_test/7);
test_vals = [horas_test dias_test semanas_test];
data_ini_test = [dados(train_vals(1)+24,1) dados(train_vals(1)+24,2) dados(train_vals(1)+24,3)];
data_fim_test = [dados(end,1) dados(end,2) dados(end,3)];
clear horas_test dias_test semanas_test

fprintf('BANCO DE DADOS (%d/%d/%d)-(%d/%d/%d)\n', data0,dataend);
fprintf('\nHoras = %d\nDias = %d\nSemanas = %d\n\n',...
    data_vals(1), data_vals(2), data_vals(3));
% Treino (70%)
fprintf('TREINAMENTO(%d/%d/%d)-(%d/%d/%d)\n', data_ini_train,data_fim_train);
fprintf('\nHoras = %d\nDias = %d\nSemanas = %d\n\n',...
    train_vals(1), train_vals(2), train_vals(3));
% Teste (30%)
fprintf('TESTE (%d/%d/%d)-(%d/%d/%d)\n', data_ini_test,data_fim_test);
fprintf('\nHoras = %d\nDias = %d\nSemanas = %d\n\n',...
    test_vals(1), test_vals(2), test_vals(3));
%% ETAPA 02 - Dados de treinamento e Estratégia de rede
% 16/04/2016 - 07/04/2017
% Y	 M	D	H	T	P	W
% Ajuste das entradas
da = str2num(teste([3])); % semanas anteriores
dp = str2num(teste([4])); % semanas posteriores
% Começo do treinamento

data0 = [16 04 2016]; 
posicao0 = find( dados(1:end,1) == data0(1) & ... % DIA
                 dados(1:end,2) == data0(2) & ... % MES
                 dados(1:end,3) == data0(3) & ... % ANO
                 dados(1:end,4) == 0);            % HORA

% N° de treinamentos
n_trains = train_vals(3) - da + 1;
% N° de horas de um treinamento
horas_trainamento = da*24;
% N° de intervalo de Semanas
n_inter = 3*24; % 3 dias, 24 horas

% vetor de posições iniciais de cada hora [TEMP PRESS WIND]
n_inf = 1;
pos_aux = 1:n_inf:n_inf*horas_trainamento;
% posição inicial e final do treinamento
ini = posicao0;
fim = horas_trainamento;
% inicialização da entrada


% Laço de criação da entrada
Treino = FORM_entradas( n_trains, ini, fim, dados, n_inf, sp, n_inter );
save('Treino.mat', 'Treino')

E0 = Treino.E;
S0 = Treino.S;
%% ETAPA 03 - Criação e treinamento da rede

if criar_rede == 1
    nnco = str2num(teste([6 7 8 9]));       % Número de neurônios da camada oculta
    net = feedforwardnet(nnco, ft);         % Definição da rede
    net.layers{1}.transferFcn = 'purelin';  % Função de ativação da camada de entrada
    net.layers{2}.transferFcn = 'purelin';  % Função de ativação da camada de saída
    net.trainParam.epochs   = 2*1e3;     % NUMERO DE EPOCAS
    net.trainParam.goal     = 1*1e-6;    % ERRO?
    net.trainParam.min_grad = 1*1e-9;    % GRADIENTE DO ERRO
    net = initlay(net); % net = init(net);
    % net.trainParam.showWindow = false; % DESABILITAR NNTRAINTOOL
    n_treinos = 15;
    [net,tr] = train(net,Treino.E,Treino.S);
    save([pwd '\Resultados - INMET\' teste '.mat'],'net', 'tr');
else
    load([pwd '\Resultados - INMET\' teste '.mat']);
end

fprintf('ETAPA 03 (OK!) - Criação e treinamento da rede \n')
%% ETAPA 04 - Teste da Rede

% posiçã inicial dos testes
posicao0_teste = train_vals(1);
% N° de testes
n_tests = 10; %test_vals(3);
% N° de horas da entrada de um teste
horas_teste_in = horas_trainamento;
horas_teste_out = sp*7*24;

% vetor de posições iniciais de cada hora [TEMP PRESS WIND]
pos_aux = 1:3:4*horas_teste_in;
% posição inicial e final do treinamento
ini = posicao0_teste;
fim = posicao0_teste + horas_teste_in-1;

% Laço de criação da entrada
Testes = FORM_entradas( n_tests, ini, fim, dados, n_inf, sp, n_inter );


for k = 1:4 %size(Testes.E,1)
    Testes.s(:,k) = sim(net,Testes.E(k,1));
end


save('Testes.mat', 'Testes')
E1 = Testes.E;
S1 = Testes.S;
n_plot = 4;

for k = 1:4
    subplot(n_plot/2,2,k);
    v_real = Testes.E(:,k); % coluna vento (17)
    v_prev = Testes.s(:,k);                % Vento previsto
    
    % Aplicação dos critérios (NIAE correlação erro)
    %[erro_p(:,k), NIAE_val(k,:) ] = NIAE_function([ha 23], h(1,ha+1:end), v_real(1,ha+1:end), v_prev);
    
    % Plotagem
    plot(v_real, 'b-o', 'Linewidth', 2); grid on; hold on
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
