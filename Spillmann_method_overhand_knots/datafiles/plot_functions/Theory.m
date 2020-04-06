

Pbar = linspace(0.01,0.2499999999,1000);
lambdaTh = Pbar;
lTh = Pbar;
Theta = linspace(0,pi,100);
for kk=1:length(Pbar)
    arg = 1./sqrt(1 - 2*Pbar(kk)*(1+cos(Theta)));
    lambdaTh(kk) = 2.0*trapz(Theta, arg);
    arg = cos(Theta)./sqrt(1 - 2*Pbar(kk)*(1+cos(Theta)));
    lTh(kk) = 2.0*trapz(Theta, arg);
end

eTh = lambdaTh + lTh;

Xth = 1/(8*sqrt(3)*pi^2) * lTh.^2./eTh;
Yth = 0.80/(96*sqrt(3)*pi^2) * lTh.^3;  %coefficient of friction

plot(Xth, Yth, 'k--');
