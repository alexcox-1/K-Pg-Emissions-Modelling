# test da gaussian perturbation
co2vals = zeros(300);
co2vals[1:75] .= 0.045;
co2vals[76:125] .= 0.005;
co2vals[106:150] .= 0.065;
co2vals[151:200] .= 0.0;
co2vals[201:250] .= 0.075;
co2vals[251:300] .= 0.035;
logco2valstest = log.(co2vals);
randsigma = rand()*length(co2valstest)/100
randsigma2 = rand()*length(co2valstest)/100
randsigma3 = rand()*length(co2valstest)/100
randmu = rand()*length(co2valstest)
randmu2 = rand()*length(co2valstest)
randmu3 = rand()*length(co2valstest)
randamplitude = randn()/normpdf(randmu, randsigma, randmu)/10
randamplitude2 = randn()/normpdf(randmu2, randsigma2, randmu2)/10
randamplitude3 = randn()/normpdf(randmu3, randsigma3, randmu3)/10
for j=1:length(co2valstest)
    logco2valstest[j] += randamplitude * normpdf(randmu, randsigma, j)
    logco2valstest[j] += randamplitude2 * normpdf(randmu2, randsigma2, j)
    logco2valstest[j] += randamplitude3 * normpdf(randmu3, randsigma3, j)
end
plot(exp.(logco2valstest))
