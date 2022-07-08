function [ erro_p, NIAE_val ] = NIAE_function(faixa_tempo, X, Y_ref, Y_comp)
%   Fun��o que calcula a Integral normalizada do erro absoluto (do ingl�s, 
%   Normalized Integral of Absolute Error ((NIAE)))

%   Defini��es
%   faixa_tempo : Vetor de duas posi��es com os valores de tempo, inicial e final;
%   X           : Vetor de tempo com passo fixo;
%   Y_ref       : Vetor da vari�vel de refer�ncia;
%   Y_comp      : Vetor da vari�vel comparada;

%Corpo da fun��o
% Defini��o do erro (passo de aquisi��o)
passo = X(2)-X(1);

% Amostras relacionadas � faixa de tempo disponibilizada
for k = 1:2
    aux = find(X < faixa_tempo(k)+passo & X > faixa_tempo(k)-passo);
    faixa_tempo_n(k) = aux(end);
end

% Defini��o do vetor de amostras
vet = faixa_tempo_n(1):faixa_tempo_n(2);

% C�lculo da area de referencia
A_ref = trapz(X(vet), Y_ref(vet));

% Calculo do m�dulo da diferen�a
diff = abs(Y_ref(vet) - Y_comp(vet));

% C�lculo da varia��o de �rea
A_delta = trapz(X(vet), diff);

% C�lculo da Integral normalizada do erro absoluto
NIAE_val = 1 - A_delta(1)/A_ref(1);

%% Corela��o
corrM = corrcoef(Y_ref, Y_comp);
corr = corrM(1,2);

%% Rela��o sinal/erro
SNR = 20*log10( rms(Y_comp(vet)) / rms(Y_ref(vet)-Y_comp(vet)) );

%% ISE Integral do erro quadratico
erro = Y_ref(vet) - Y_comp(vet);
erro2 = Y_comp(vet) - (circshift(Y_comp(vet)',-1))';

erro_p = (abs(erro)./Y_ref(vet))';
erro_p(~isfinite(erro_p)) = 0;
% ISE = trapz(X(vet), erro.^2);

%% Integral do erro absoluto
%IAE = trapz(X(vet), abs(erro));

%% EMQ erro m�dio quadratico
NMSE = (1/length(erro))*sum(erro.^2)*(100/rms(Y_ref));

%% MAPE (�Mean Absolute Percentage Error�); 
MAPE = 100/length(erro_p)*mean(erro_p);

%%  RMSE (�Root Mean Square Error�)
RMSE = sqrt(sum(erro.^2)/length(erro));

%%  O coeficiente U de Theil. 
U = sqrt(sum(erro.^2)) / sqrt(sum((  erro2  ).^2));

%% Composi��o do vetor respostas
% corr ISE IAE
NIAE_val = [MAPE RMSE U NIAE_val ];%corr NMSE ];

% NIAE_val(end+1, :) = [mean(NIAE_val(end,:))];

%% xcorr
% [r, lags] = xcorr(Y_ref(vet), Y_comp(vet), 'coeff');
% NIAE_val = [NIAE_val max(r)];

%% corr2
% NIAE_val = [NIAE_val corr2(Y_ref(vet(1:end-10)), Y_comp(vet(1:end-10)))];

end

