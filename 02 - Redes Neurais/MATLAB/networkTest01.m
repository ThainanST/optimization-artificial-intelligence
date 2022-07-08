function [ S ] = networkTest01( net, E )

num = length(E);
S = zeros(1,num);
for k = 1:num
    S(k) = sim(net,E(:,k));
end

end

