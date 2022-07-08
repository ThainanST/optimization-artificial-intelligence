% Exemplo de aplicação de Redes Neurais Artificiais

clear all
clc

% Conjunto de treinamento (tr)
E1 = [-0.5 -0.375 -0.25 -0.125 0 0.125 0.25 0.375 0.5]; % X
S1 = [-1.3268E-6 -0.35355 -0.5 -0.35355 0 0.35355 0.5 0.35355 1.3268E-6]; % Y


% Conjunto de teste (te)
E2 = [-0.47747 -0.381972 -0.254648 -0.111408 -0.015916 0.111409 0.238733 0.33342256 0.5252117]; % X
S2 = [-1.E-6 -0.4 -0.56 -0.4 0 0.3 0.6 0.4 1.2E-6]; % Y


ft='trainrp';

nnco = 3; % Número de neurônios da camada oculta
nncs = 1; % Número de neurônios da camada de saída

fs  = 'tansig';
fco = 'tansig';

[E1N,mine,maxe,S1N,mins,maxs] = premnmx(E1,S1);

E2N = tramnmx(E2,mine,maxe);
S2N = tramnmx(S2,mins,maxs);

net = newff(minmax(E1N),[nnco,nncs],{fco,fs},ft); % 1 camada oculta !

test.P = E2N;
test.T = S2N;

net.trainParam.show   = 100;
net.trainParam.epochs = 100000;
net.trainParam.goal   = 0.005;

net=init(net);

[net,tr]=train(net,E1N,S1N,[],[]);
net.IW{1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(tr.epoch, tr.perf);
legend('Treinamento',-1);
ylabel('erro quadratico'); xlabel('Epoca');

R1N = sim(net,E1N); % Saída da rede para o conjunto de treinamento

[R1] = postmnmx(R1N,mins,maxs);

R2N = sim(net,E2N); % Saída da rede para o conjunto de teste

[R2]=postmnmx(R2N,mins,maxs);

plot(E2,R2);
