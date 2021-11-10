function saveFigasPDF(locs,name)
for i = 1:length(locs)
    exportgraphics(gcf,strcat(locs{i}, name, '.pdf'),...
        'ContentType','vector');
end
end