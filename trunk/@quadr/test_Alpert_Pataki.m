%
% Test the QuadNodeInterval function
%

IN = {
% smooth functions
    @(x) exp(x),                  0.0,  1.0,  1,1, 1.7182818284590452354,
    @(x) exp(-x.^2),              0.0, 10.0,  1,0, 0.88622692545275801365,
    @(x) exp(-x.^2),            -10.0, 10.0,  0,0, 2*0.88622692545275801365,
% square root singularities
    @(x) 1./sqrt(x),              0.0,  1.0,  2,1, 2.0,
    @(x) 1./sqrt(x.*(1-x)),       0.0,  1.0,  2,2, 3.1415926535897932385,
    @(x) exp(-x.^2)./sqrt(x),     0.0, 10.0,  2,1, 1.8128049541109541559,
% log singularities
    @(x) log(x)			  0.0,  1.0,  3,1, -1.0,
    @(x) exp(x).*log(x)		  0.0,  1.0,  3,1, -1.3179021514544038949,
    @(x) exp(-x.^2).*log(x)	  0.0, 10.0,  3,1, -.87005772672831550673,
% mixed ones
    @(x) log(x)./sqrt(1-x)        0.0,  1.0,  3,2, -1.2274112777602187623,
    };


N = 40;
order = 10;
Kmax = size(IN,1);

for K=1:Kmax

    fn = IN{K,1};
    a = IN{K,2};
    b = IN{K,3};
    corra = IN{K,4};
    corrb = IN{K,5};
    Iexact = IN{K,6};
    
    [Ax, Aw] = QuadNodesInterval(a, b, N, 0, corra, corrb, order)
    I = sum(fn(Ax).*Aw);
    Ierr = I - Iexact;
    fprintf('Case:%02d,  I=%23.16e,  Iexact=%23.16e,  Ierr=%13.6e\n', K, I, Iexact, Ierr);
    
end
