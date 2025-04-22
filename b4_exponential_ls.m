function ls = b4_exponential_ls(params,x,y)

a = params(1);
b = params(2);
r = params(3);

predy = b5_exponential(x,a,b,r);
ls = sum( (y-predy).^2);
