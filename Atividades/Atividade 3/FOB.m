function [ fob ] = FOB( Dados_prob, MORCEGO )

y1 = Dados_prob.x1x2y([1 2 3],:)'*MORCEGO(1,1:Dados_prob.d)';
erro = y1 - Dados_prob.x1x2y(end,:)';
fob = sqrt(erro.^2)/length(erro);
fob = mean(fob);

end