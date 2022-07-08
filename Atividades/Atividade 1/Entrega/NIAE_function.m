function [ erro_p, NIAE_val ] = NIAE_function(faixa_tempo, X, Y_ref, Y_comp)
%   Função que calcula a Integral normalizada do erro absoluto (do inglês, 
%   Normalized Integral of Absolute Error ((NIAE)))

%   Definições
%   faixa_tempo : Vetor de duas posições com os valores de tempo, inicial e final;
%   X           : Vetor de tempo com passo fixo;
%   Y_ref       : Vetor da variável de referência;
%   Y_comp      : Vetor da variável comparada;

%Corpo da função
% Definição do erro (passo de aquisição)
passo = X(2)-X(1);

% Amostras relacionadas à faixa de tempo disponibilizada
for k = 1:2
    aux = find(X < faixa_tempo(k)+passo & X > faixa_tempo(k)-passo);
    faixa_tempo_n(k) = aux(end);
end

% Definição do vetor de amostras
vet = faixa_tempo_n(1):faixa_tempo_n(2);

% Cálculo da area de referencia
A_ref = trapz(X(vet), Y_ref(vet));

% Calculo do módulo da diferença
diff = abs(Y_ref(vet) - Y_comp(vet));

% Cálculo da variação de área
A_delta = trapz(X(vet), diff);

% Cálculo da Integral normalizada do erro absoluto
NIAE_val = 1 - A_delta(1)/A_ref(1);

%% Corelação
corrM = corrcoef(Y_ref, Y_comp);
corr = corrM(1,2);

%% Relação sinal/erro
SNR = 20*log10( rms(Y_comp(vet)) / rms(Y_ref(vet)-Y_comp(vet)) );

%% ISE Integral do erro quadratico
erro = Y_ref(vet) - Y_comp(vet);
erro2 = Y_comp(vet) - (circshift(Y_comp(vet)',-1))';

erro_p = (abs(erro)./Y_ref(vet))';
erro_p(~isfinite(erro_p)) = 0;
% ISE = trapz(X(vet), erro.^2);

%% Integral do erro absoluto
%IAE = trapz(X(vet), abs(erro));

%% EMQ erro médio quadratico
NMSE = (1/length(erro))*sum(erro.^2)*(100/rms(Y_ref));

%% MAPE (“Mean Absolute Percentage Error”); 
MAPE = 100/length(erro_p)*mean(erro_p);

%%  RMSE (“Root Mean Square Error”)
RMSE = sqrt(sum(erro.^2)/length(erro));

%%  O coeficiente U de Theil. 
U = sqrt(sum(erro.^2)) / sqrt(sum((  erro2  ).^2));

%% Composição do vetor respostas
% corr ISE IAE
NIAE_val = [MAPE RMSE U NIAE_val ];%corr NMSE ];

% NIAE_val(end+1, :) = [mean(NIAE_val(end,:))];

%% xcorr
% [r, lags] = xcorr(Y_ref(vet), Y_comp(vet), 'coeff');
% NIAE_val = [NIAE_val max(r)];

%% corr2
% NIAE_val = [NIAE_val corr2(Y_ref(vet(1:end-10)), Y_comp(vet(1:end-10)))];

end

