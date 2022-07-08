%% Rede Hopfield
clc
clf
close all
clear all

%% Inicialização de parametros
Convergiu = false;	    			% indica se a rede está treinada
FimTreinamento = false;		        % indica fim do treinamento
ErroEnergia=1e-5;                   % energia do erro entre o padrao de saida na época atual e anterior
MaxEpocas=10000;                    % quantidade maxima de épocas de treinamento
Epc=0;                              % indica a época corrente. Inicia na época 0
P=4;                                % qtd de Padroes a serem armazenados na matriz de pesos
linha=9;                            % qtd linha em cada imagem
coluna=5;                           % qtd colunas em cada imagem
Neuronios=45;                       % qtd. de elementos dos vetores de entrada e tambem a quantidade de neuronios da rede
M=4;                                % m imagens corrompidas com ruido para cada padrao armazenado.
%%
nivelRuido=20/100;
QtdAlterar=floor(Neuronios*nivelRuido);
PixelAlterar=0;
ImgRecuperada=zeros(linha,coluna);
Separador=ones(size(ImgRecuperada));
Imagem=cell(linha,coluna,P);

%% imagem do 1
Imagem(1)={[       -1 -1 -1 -1 -1
                   -1 -1 -1 -1 -1
                    0  0 -1  0  0
                    0  0 -1  0  0
                    0  0 -1  0  0
                    0  0 -1  0  0
                    0  0 -1  0  0
                   -1 -1 -1 -1 -1
                   -1 -1 -1 -1 -1  ]};
% imagem do A
Imagem(2)={[       -1 -1 -1 -1 -1
                   -1  0  0  0 -1
                   -1  0  0  0 -1
                   -1  0  0  0 -1
                   -1 -1 -1 -1 -1
                   -1 -1 -1 -1 -1
                   -1  0  0  0 -1
                   -1  0  0  0 -1
                   -1  0  0  0 -1  ]};
% imagem do 3
Imagem(3)={[       -1 -1 -1 -1 -1
                   -1 -1 -1 -1 -1
                    0  0  0 -1 -1
                    0  0  0 -1 -1
                   -1 -1 -1 -1 -1
                    0  0  0 -1 -1
                    0  0  0 -1 -1
                   -1 -1 -1 -1 -1
                   -1 -1 -1 -1 -1]};
% imagem do 4
Imagem(4)={[       -1 -1  0 -1 -1
                   -1 -1  0 -1 -1
                   -1 -1  0 -1 -1
                   -1 -1  0 -1 -1
                   -1 -1 -1 -1 -1
                   -1 -1 -1 -1 -1
                    0  0  0 -1 -1
                    0  0  0 -1 -1
                    0  0  0 -1 -1]};
                
for p=1:P
    for l = 1:9
        for c = 1:5
            if Imagem{p}(l,c) == 0
                Imagem{p}(l,c) = 1;
            end
        end
    end
end

%% gerar matriz de pesos
W=zeros(Neuronios,Neuronios);
for p=1:P
    W=W+Imagem{p}(:)*Imagem{p}(:)';% - eye(Neuronios,Neuronios);
end
 W=W/Neuronios;

%% gerar imagens com ruido
ImagemRuido=Imagem;
% %% gera Imagens com ruido
for p=1:P
    for n=1:QtdAlterar
        PixelAlterar=floor(Neuronios*(rand))+1;
        ImagemRuido{p}(PixelAlterar)=-ImagemRuido{p}(PixelAlterar);
    end
    qtdPixeisAlterados=(sum(abs( (Imagem{p}(:) ~= ImagemRuido{p}(:) ) )));
end

%% 
for p=1:P
    Epc=0;                                % inicializa epoca para cada imagem com ruido apresentado
    Convergiu = false;	    		      % inicializa condição sinalizadora de que a rede está treinada
    FimTreinamento = false;		      
    v_atual=ImagemRuido{p};
    while ~(Convergiu || FimTreinamento)
        Epc = Epc + 1;
        v_anterior=v_atual;
        u=W*v_anterior(:);
        v_atual=sign(u);
        dist=0.5*( norm(v_anterior(:) - v_atual(:)))^2;        
        Convergiu             = ( dist  < ErroEnergia    );
        FimTreinamento = ( Epc > MaxEpocas );
    end
    
    % recompoe a imagem como matriz
    for c=1:coluna
        for  l=1:linha
            ImgRecuperada(l,c)=v_atual(l+(c-1)*linha);
        end
    end

    figure('name',['Resultados ' num2str(p) ],'NumberTitle','off');
    imshow( [ Imagem{p}        Separador ...
              ImagemRuido{p}  Separador ...
              ImgRecuperada]  ,'InitialMagnification','fit');
          drawnow
end