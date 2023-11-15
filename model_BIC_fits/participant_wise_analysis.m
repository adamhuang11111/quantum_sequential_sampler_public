nLLS_diff = zeros(1162,1);
BIC_diff = zeros(1162,1);
chi2_test = zeros(1162,1);
for i = 1:1162
    nLL_noint = nLLS_noint(i);
    nLL_int = nLLS_int(i);
    nLLS_diff(i) = nLL_noint - nLL_int;
    BIC_diff(i) = nLL_noint - nLL_int - log(78);
    chi2_test(i) = chi2pdf(nLLS_diff(i),1);
end