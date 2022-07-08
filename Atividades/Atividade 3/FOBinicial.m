function [ Solucao ] = FOBinicial( Dados_prob, Dados_simu, Solucao )

for k = 1:Dados_simu.N
    y1 = Dados_prob.x1x2y([1 2 3],:)'*Solucao.MORCEGOS(k,1:Dados_prob.d)';
    erro = y1 - Dados_prob.x1x2y(end,:)';
    fob = sqrt(erro.^2)/length(erro);
    fob = mean(fob);
end


end

