function mp=norminv(p,mu,sigma)
mp=sigma.*(-sqrt(2)*erfcinv(2*p))+mu;
return