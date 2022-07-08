% Exemplo de aplicação de Redes Neurais Artificiais
% RAIZ QUADRADA DE UM NÚMERO

clear all
close all

clc

%% CONJUNTOS DE TREINAMENTO
% Conjunto de treinamento (tr)
% AS ENTRADAS DEVEM ESTAR EM COLUNA
E1 = [ 0 1 2 3 4 5 6 7 8 9 ];
% S1 = [0 1 1.4142 1.7321 2 2.2361 2.4495 2.6458 2.8284 3];
S1 = sqrt(E1);

V1 = [E1' S1']

% ESPECIFICAÇÃO DE TREINO
ft = 'trainrp';
ft = 'trainlm';
% train: TREINO
% r    : RESILIENCIA (??)
% p    : BACKPROPAGATION

%% CRIAÇÃO DA REDE

% NUMERO DE NEURÔNIOS
nnco = 5; % Número de neurônios da camada oculta
nncs = 1; % Número de neurônios da camada de saída

% FUNÇÃO DE ATIVAÇÃO
% fs  = 'tansig';
% fco = 'tansig';
fs  = 'purelin';
fco = 'purelin';

% CRIAÇÃO DA REDE
net = newff(minmax(E1),[nnco,nncs],{fco,fs},ft); % 1 camada oculta !

%% VISUALIZAÇÃO DE PARAMETROS DA REDE
net.layers{1};
net.IW{1};
net.LW{2};
net.b{1};
net.b{2};

%% TESTE DA REDE
S2 = sim(net,E1);
% S2 = round(S2);
V2 = [E1' S2']
    
%% TREINAMENTO DA REDE

net.trainParam.epochs = 1000;    % NUMERO DE EPOCAS
net.trainParam.goal   = 0.01;    % ERRO?
net.trainParam.min_grad = 1e-5;  % GRADIENTE DO ERRO
net = init(net);

[net,tr]=train(net,E1,S1,[],[]);

%% VISUALIZAÇÃO DOS RESULTADOS
% plot(tr.epoch, tr.perf);
% ylabel('erro quadratico'); xlabel('Epoca');

S3 = sim(net,E1);
% S3 = round(S3);
V3 = [E1' S3']
    