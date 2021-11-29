function aerofoil = polars()
% extract aerofoil input data

S = dir('..\polars\*.txt');
S = S(1:6,1); % remove 'hints.txt'
aerofoil = cell(2,length(S));
aerofoil(1,:) = {S.name};
vars = {'alpha','cl','cd','fsst','clinv','clfs','cm'};

for i = 1:size(aerofoil,2)
    aerofoil{2,i} = readtable("..\polars\" + aerofoil{1,i});
    aerofoil{2,i}.Properties.VariableNames = vars;
    aerofoil{1,i} = erase(aerofoil{1,i},"_ds.txt");
end
end