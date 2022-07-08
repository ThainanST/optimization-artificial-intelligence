% 23/04/2017
% Técnicas Inteligentes - Bat Algorithm
% Prof. Ivo Chaves
% Análise de todas as possibilidades

clear all
close all
clc
N = 21;

Dados_prob.A = linspace(-10,10,N);
Dados_prob.B = Dados_prob.A;
Dados_prob.C = Dados_prob.A;

% Dados_prob.A = linspace(3.6316,3.8421,N);   % linspace(3.2,3.5,N);
% Dados_prob.B = linspace(3.4211,3.6316,N);   % linspace(3.4,3.6,N);
% Dados_prob.C = linspace(-6.5789,-6.3684,N); % linspace(-5,7,N);

% Dados_prob.A = linspace(3.6981,3.7202,N);   % linspace(3.2,3.5,N);
% Dados_prob.B = linspace(3.4543,3.4765,N);   % linspace(3.4,3.6,N);
% Dados_prob.C = linspace(-6.3906,-6.3684,N); % linspace(-5,7,N);


[Dados_prob.A' Dados_prob.B' Dados_prob.C'] 

fprintf('*********************\n');
fprintf('Solução do problema, %d possibilidades\n', length(Dados_prob.A)^3);
fprintf('*********************\n');

Dados_prob.x1x2y = [1  1 1 1 1 1 1 1; % A
                   -1  0 1 2 4 5 5 6; % B
                   -2 -1 0 1 1 2 3 4; % C
                    13 11 9 4 11 9 1 -1];
Dados_prob.d = 3;                

h = waitbar(0,'Criando as possibilidades... (0)');
n = 1;
for r = 1:length(Dados_prob.A)
    for s = 1:length(Dados_prob.B)
        for t = 1:length(Dados_prob.C)
            Dados_prob.ABC(n, [1 2 3]) = [Dados_prob.A(r) Dados_prob.B(s) Dados_prob.C(t)];
            n = n + 1;
        end
    end
    waitbar(r/length(Dados_prob.A),h);
    ntext = ['Criando as possibilidades...' num2str(100*(r)/length(Dados_prob.A), '(%0.f)')];
    set( get(findobj(h,'type','axes'),'title'), 'string', ntext);
end, clear r s t n

close(h)

h = waitbar(0,'Testando FOB... (0)');

for m = 1:length(Dados_prob.ABC)
    Dados_prob.ABC(m,Dados_prob.d+1) = FOB( Dados_prob, Dados_prob.ABC(m,1:Dados_prob.d) );
    
    waitbar(m/length(Dados_prob.ABC),h);
    ntext = ['Testando FOB...' num2str(100*(m)/length(Dados_prob.ABC), '(%0.f)')];
    set( get(findobj(h,'type','axes'),'title'), 'string', ntext);
end, clear m

close(h)

% [r, index] = sort(Solucao.Teste(:,end));
% Solucao.Teste_ord = Solucao.Teste(index,:);
% Dados_prob.ABC_ord = Dados_prob.ABC(index,:);

[index, r2] = find(min(Dados_prob.ABC(:,end)) == Dados_prob.ABC(:,end));


fprintf('A = %.4f\nB = %.4f\nC = %.4f\n',Dados_prob.ABC(index(1),1:Dados_prob.d))
fprintf('*********************\n');
fprintf('Menor erro médio = %.4f\n', Dados_prob.ABC(index(1),end));
fprintf('*********************\n');

