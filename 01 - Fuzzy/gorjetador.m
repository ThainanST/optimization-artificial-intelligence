%===============================================================================%
%                                                                               %    
%                       UNIVERSIDADE FEDERAL DE JUIZ DE FORA - UFJF             %
%                          PÓS-GRADUAÇÃO EM ENGENHARIA ELETRICA                 %
%                                 TÉCNICAS DE OTIMIZAÇÃO                        %
%                                 SISTEMAS INTELIGENTES                         %
%                                                                               %
% Professor : Leonardo Willer de Oliveira                                       %
%                                                                   20/03/2017  %
%===============================================================================%


clear all
clc


dados=[6 8; 1 5];


%-------------------------------------------
%  Execucao do Sistema de Inferencia FUZZY 
%-------------------------------------------

a=readfis('gorjeta.fis');

Resultado = evalfis(dados,a)
