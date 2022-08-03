function output = levy()
    beta = 1.5;    
    k1 = 1;
    k2 = 1.849902656990532;
    
    x = (rand(1)*2-1)*2;
    output = normpdf(x,k2)/(abs(randn(1)))^beta;
end