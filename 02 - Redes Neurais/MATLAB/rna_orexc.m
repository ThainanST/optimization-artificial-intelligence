% Exemplo de aplicação de Redes Neurais Artificiais

clear all
clc

%% CONJUNTOS DE TREINAMENTO
% Conjunto de treinamento (tr)
% AS ENTRADAS DEVEM ESTAR EM COLUNA
E1 = [ 0 0 1 1
       0 1 0 1];

%   RESPECTIVAS SAÍDAS (PORTA OU)
S1 = [0 1 1 1];


% ESPECIFICAÇÃO DE TREINO
ft = 'trainrp';
% train: TREINO
% r    : RESILIENCIA (??)
% p    : BACKPROPAGATION

%% CRIAÇÃO DA REDE

% NUMERO DE NEURÔNIOS
nnco = 1; % Número de neurônios da camada oculta
nncs = 1; % Número de neurônios da camada de saída

% FUNÇÃO DE ATIVAÇÃO
% fs  = 'tansig';
% fco = 'tansig';
fs  = 'purelin';
fco = 'purelin';

% CRIAÇÃO DA REDE
net = newff(minmax(E1),[nnco,nncs],{fco,fs},ft);
% Syntax
%      net = newff(P,T,S)
%      net = newff(P,T,S,TF,BTF,BLF,PF,IPF,OPF,DDF)
%    Description
%      newff(P,T,S) takes,
%        P  - RxQ1 matrix of Q1 representative R-element input vectors.
%        T  - SNxQ2 matrix of Q2 representative SN-element target vectors.
%        Si  - Sizes of N-1 hidden layers, S1 to S(N-1), default = [].
%              (Output layer size SN is determined from T.)
%      and returns an N layer feed-forward backprop network.
%  
%      newff(P,T,S,TF,BTF,BLF,PF,IPF,OPF,DDF) takes optional inputs,
%        TFi - Transfer function of ith layer. Default is 'tansig' for
%              hidden layers, and 'purelin' for output layer.
%        BTF - Backprop network training function, default = 'trainlm'.
%        BLF - Backprop weight/bias learning function, default = 'learngdm'.
%        PF  - Performance function, default = 'mse'.
%        IPF - Row cell array of input processing functions.
%              Default is {'fixunknowns','remconstantrows','mapminmax'}.
%        OPF - Row cell array of output processing functions.
%              Default is {'remconstantrows','mapminmax'}.
%        DDF - Data division function, default = 'dividerand';
%      and returns an N layer feed-forward backprop network.
%load net

%% VISUALIZAÇÃO DE PARAMETROS DA REDE

net.layers{1}; % CAMADAS DA REDE
net.IW{1};     % GANHOS DE ENTRADA
net.LW{2};     % GANHO INTERMEDIÁRIO
net.b{1} ;     % BIAS DO NEURÔNIO DE ENTRADA
net.b{2};      % BIAS DO NEURÔNIO DE SAÍDA

%% TESTES NA REDE

S = networkTest01( net, E1 );
S = round(S);
V = [E1; S]


%% TREINAMENTO DA REDE

% DEFINIÇÃO DE PARAMETROS
net.trainParam.epochs = 1000;    % NUMERO DE EPOCAS
net.trainParam.goal   = 0.01;    % ERRO?
net.trainParam.min_grad = 1e-5;  % GRADIENTE DO ERRO
% INICALIZAÇÃO DA REDE
net = init(net);

% TREINO DA REDE
[net,tr] = train(net,E1,S1,[],[]);

%% VISUALIZAÇÃO DOS DADOS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(tr.epoch, tr.perf, tr.epoch, tr.vperf);
% plot(tr.epoch, tr.vperf);
plot(tr.epoch, tr.perf);
% hold on
% plot(tr.epoch, tr.vperf);
ylabel('erro quadratico'); xlabel('Epoca');

S = networkTest01( net, E1 );
S = round(S);
V = [E1; S]

