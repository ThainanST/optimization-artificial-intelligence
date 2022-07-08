% Atividade Técnicas Inteligentes - RNA
% Prof. Leonardo Willer de Oliveira

% Thainan Santos Theodoro
% Juiz de Fora, 27 de Março de 2017
% https://www.mathworks.com/help/nnet/ref/feedforwardnet.html

% ESTRATÉGIA: 
% 3 (três) horas anteriores determinam a 1 (primeira) hora
%  seguinte.
% Usando feedforwardnet
% Banco de dados de itapipoca 
% http://www.inmet.gov.br/portal/index.php?r=estacoes/estacoesAutomaticas

clear all
close all
clc
teste = 'JF21-05'; % A
% teste = 'JF21-10'; % B
% teste = 'JF21-15'; % C
% teste = 'JF21-20'; % D
% teste = 'JF22-05'; % E
% teste = 'JF22-10'; % F
% teste = 'JF22-15'; % G
% teste = 'JF22-20'; % H
% teste = 'JF31-05'; % I
% teste = 'JF31-10'; % J
% teste = 'JF31-15'; % K
% teste = 'JF31-20'; % L

% teste = 'JF32-05'; % M
% teste = 'JF32-10'; % N
% teste = 'JF32-15'; % O
% teste = 'JF32-20'; % P
% teste = 'JF41-05'; % Q
% teste = 'JF41-10'; % R
% teste = 'JF41-15'; % S
% teste = 'JF41-20'; % T
% teste = 'JF42-05'; % U
% teste = 'JF42-10'; % V
% teste = 'JF42-15'; % X
% teste = 'JF42-20'; % W

criar_rede = 1;
% ft = 'trainrp';
ft = 'trainlm';
% ft = 'trainscg';
% ft = 'trainbr';


fprintf('Teste: %s \nCriar rede = %d \nTreino = %s\n\n', teste,  criar_rede, ft)

%% ETAPA 01 - Leitura dos dados
dados = xlsread('\Dados\juizdefora');
% dados = xlsread('\Dados\itapipoca');
fprintf('ETAPA 01 (OK!) - Leitura dos dados \n')

%% ETAPA 02 - Dados de treinamento e Estratégia de rede


ha = str2num(teste([3])); % horas anteriores
hp = str2num(teste([4])); % horas posteriores
% data0 = [30 03 2016]; % [ano mes dia] Dia do treinamento
data0 = [2 04 2016]; 
posicao0 = find( dados(1:end,1) == data0(3) & ... % ANO
                 dados(1:end,2) == data0(2) & ... % MES
                 dados(1:end,3) == data0(1) )';   % DIA
% Configuração da entrada
horas_train = 24*7*4*8; % 24 672
for k = 1:horas_train-ha % 24 horas
    ini = posicao0(1)-1+k;
    fim = ini + ha - 1;
    E0(:,k) = dados(ini:fim,17); % coluna vento (17)
end, clear k ini fim

% Configuração da saída
for k = 1:length(E0)
    ini = posicao0(1)+ha+hp-2+k;
    fim = ini + hp - 1;
    S0(:,k) = dados(ini:fim,17); % coluna vento (17)
end, clear k ini fim
% V0 = [E0; S0]';

fprintf('ETAPA 02 (OK!) - Dados de treinamento \n')

%% ETAPA 03 - Criação e treinamento da rede

if criar_rede == 1
    nnco = str2num(teste([6 7]));       % Número de neurônios da camada oculta
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

    net = feedforwardnet(nnco, ft); % Definição da rede
    
    % Função de ativação da camada de entrada
    %      tansig Symmetric sigmoid transfer function
    %      purelin Linear transfer function
    % 'tansig' 'purelin'
    net.layers{1}.transferFcn = 'purelin'; 
    net.layers{2}.transferFcn = 'purelin';
    
    % net.performFcn = 'mse';  %Name of a network performance function %type help nnperformance
    % view(net)
%     net = init(net)
%     net = initzero(net, [1 1]); % initzero(S,PR)
    
    
    % parâmetros da rede 
    net.trainParam.epochs   = 1*1e3;     % NUMERO DE EPOCAS
    net.trainParam.goal     = 1*1e-6;    % ERRO?
    net.trainParam.min_grad = 1*1e-9;    % GRADIENTE DO ERRO
    
    % treino da rede
    % net.trainParam.showWindow = false; % DESABILITAR NNTRAINTOOL
%     net = init(net);
    net = initlay(net);
    [net,tr] = train(net,E0,S0);
    %save([pwd '\Resultados - INMET\' teste '.mat'],'net', 'tr');
else
    %load([pwd '\Resultados - INMET\' teste '.mat']);
end

% Plotagem de gráfico de performance do treinamento
if 1 == 0
    plotperf(tr); hold on
    plot(tr.best_epoch, tr.best_perf,  'ro');
    plot(tr.best_epoch*[1 1+1e-9], [1e-1 max(tr.perf)], 'g--')
    ylim([1e-1 max(tr.perf)])
    legend(['Treino'], ['Validação = ' num2str(tr.best_vperf) ], ['Teste = ' num2str(tr.best_tperf)]) 

end


% Parametros da rede
net.layers{1};
net.IW{1};
net.LW{2}';
net.b{1};
net.b{2};

fprintf('ETAPA 03 (OK!) - Criação e treinamento da rede \n')

%% ETAPA 04 - Teste da Rede

% Dias do teste (sempre par) 
% dia do teste(2 04 2016) + 8meses
datas = [ 10 01 2017
          20 01 2017
          10 02 2017
          20 02 2017];

[n_plot, aux] = size(datas); % Número de plots (sempre par)
posicao = zeros(n_plot, 24); % Inicialização

% Posição dos dados de vento nos dados
for k = 1:n_plot
    posicao(k,:) = find( dados(1:end,1) == datas(k,3) & ... % ANO
                         dados(1:end,2) == datas(k,2) & ... % MES
                         dados(1:end,3) == datas(k,1) )';   % DIA
end, clear k

% Configuração das entradas e aplicação da rede
for m = 1:n_plot
    % Configuração das entradas
    for k = 1:length(posicao)-ha
        ini = posicao(m,1)-1+k;
        fim = ini + ha - 1;
        E(:,k) = dados(ini:fim,17); % coluna vento (17)
    end
    % Aplicação da rede
    aux = sim(net,E);
    for n = 1:size(aux, 1)
        S(m,:,n) = aux(n,:);
    end
    clear k n ini fim aux
end

fprintf('ETAPA 04 (OK!) - Teste da Rede \n')

%% ETAPA 05 - visualização dos resultados

% Vetor de horas reais
h = dados(posicao0,4)'; % coluna de horas (4);

% Vetor de horas previstas
h_prev = ha:ha-1+length(S); % coluna de horas (4)

% Visualização
close all
% screenposition = get(fig,'Position');
%screenposition = [319   358   854   586]; %391   371   993   333 para 2 graficos
%fig = figure('Position', screenposition);
% fig = figure;
for k = 1:n_plot
    subplot(n_plot/2,2,k);
    v_real = dados(posicao(k,:),17)'; % coluna vento (17)
    v_prev = S(k,:,1);                % Vento previsto
    
    % Aplicação dos critérios (NIAE correlação erro)
    [erro_p(:,k), NIAE_val(k,:) ] = NIAE_function([ha 23], h(1,ha+1:end), v_real(1,ha+1:end), v_prev);
    
    % Plotagem
    plot(h     , v_real, 'b-o', 'Linewidth', 2); grid on; hold on
    plot(h_prev, v_prev, 'r-o', 'Linewidth', 2); ...
    xlim([0 23]), ylim([0 10]);
    title([ '\bf(' teste ')(' ft ')Previsão do dia ' num2str(datas(k,1)), '/', num2str(datas(k,2)),'/', num2str(datas(k,3)) ]);
%     title([ num2str(datas(k,1)), '/', num2str(datas(k,2)),'/', num2str(datas(k,3)), ...
%         '; MAPE = ', num2str(NIAE_val(k,1)), '; U = ', num2str(NIAE_val(k,3))]);
    if k == 1
    %text(0.5, 9, ['\bf' 'CONFIG: ' '(' num2str(ha) 'x' ...
     %   num2str(hp) ') (' num2str(nnco) 'n)']); % Configuração do treinamento
                                                % (ha x hp) Neuronios n
    % text(0.5, 7,['\bf' 'DIA DO TESTE: ' , num2str(data0(1)), '/', num2str(data0(2)),'/', num2str(data0(3))]);
    end
    legend('Real','Previsto');
    xlabel('Tempo (hora)');
    ylabel('Velocidade (m/s)')
    
    % Plotagem se houver mais de uma saída
    if hp > 1
        v_prev2 = S(k,:,2);
        h_prev2 = ha+hp-1:ha+length(S);
        plot(h_prev2, v_prev2, 'm-o', 'Linewidth', 2);
        legend('Real','1° hora','2° hora');
    end
    
end

% Inclusão de última linha com valores médios dos dias selecionados
NIAE_val(end+1, :) = mean(NIAE_val, 1)

% Visualição dos valores médios selecionados
% [MAPE RMSE U NIAE_val corr NMSE ]
% [ NIAE_val(end, [1 2])'; 100*NIAE_val(end, [5 6])' ];














