function [out,a] = flattenTip(var1,var2,rotor,r_R)
a = var1( find( (rotor.r/rotor.R) > r_R, 1 ) );
var1((rotor.r/rotor.R) > r_R) = ones(1,length(var2((rotor.r/rotor.R) > r_R))) .* a;
out = var1;
end