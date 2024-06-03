function new_beta = betaSpline(redesign,result,span,ratios,beta_max)

r         = redesign.r;
R         = redesign.R;
beta      = result.beta;
grad_beta = gradient(beta);

start  = find(r/R >= span(1), 1, 'first');
middle = find(r/R >= span(2), 1, 'first');
final  = find(r/R <= span(3), 1, 'last' );
pos    = [start middle final];

y_middle = beta(final) + (beta(start) - beta(final)) * ratios(1);
y        = [beta(start), y_middle, beta(final)];
beta_spl = spline(r(pos), [grad_beta(start) y grad_beta(final)*ratios(2)], r(start:final));

new_beta = beta;
new_beta(start:final) = beta_spl;

new_beta(new_beta > deg2rad(beta_max)) = deg2rad(beta_max);

end