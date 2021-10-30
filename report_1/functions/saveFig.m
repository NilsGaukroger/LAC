function saveFig(locs,name,type)
for i = 1:length(locs)
    if type == "pdf"
        exportgraphics(gcf,strcat(locs{i}, name, '.pdf'),...
            'ContentType','vector');
    elseif type == "png"
        exportgraphics(gcf,strcat(locs{i}, name, '.png'));
    end
end
end