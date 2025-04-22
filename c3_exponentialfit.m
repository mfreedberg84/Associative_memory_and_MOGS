function [bestparams, rr, minls, exitflag] = c3_exponentialfit(xdata,ydata)

%if no exceptions file is specified use this.
starting_a = max(ydata);
starting_b = max(ydata) - min(ydata);
starting_r = 1;

% boundariES on A, B and R
lb_a = 0;
ub_a = max(ydata);
    
lb_b = 0;
ub_b = 10;

lb_r = 0;
ub_r = 10;

%set up start parameters and options
basestartparams = [starting_a starting_b starting_r];
 
options =optimset('Display','off','Tolfun',1e-30,'TolX',1e-8,'MaxIter',200000,'MaxFunEvals',400000);
warning off all

bestfval = inf;

startparams = basestartparams;
[params, fval, exitflag] = fmincon(@(x)c4_exponential_ls(x,xdata,ydata),startparams,[],[],[],[],[lb_a lb_b lb_r],[ub_a ub_b ub_r],[],options);

if fval < bestfval
	bestparams = params;
	bestexitflag = exitflag;
end

predy = c5_exponential(xdata, bestparams(1), bestparams(2), bestparams(3));
minls = sum((predy - ydata).^2);
               
rr = corrcoef(predy, ydata);
rr = rr(1,2);

exitflag = bestexitflag;

                        
