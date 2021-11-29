function [files, channels, data] = load_stats(stat_file,turbine)
% Load the statistics from the text file
% Load the raw data
raw        = readtable(stat_file,'ReadVariableNames',true,'VariableNamingRule',"preserve");
channels   = cellfun(@str2double,raw.Properties.VariableNames(2:end));
n_files    = size(raw,1);
n_channels = length(channels);

wsps = NaN(1,n_files);
for i = 1:n_files
    if turbine == 1
        wsps(i) = str2double(regexp(raw.File{i},'(?<=turb_)\d+','match'));
    elseif turbine == 2
        wsps(i) = str2double(regexp(raw.File{i},'(?<=dl_)\d+','match'));
    end
end

[wsps_unique,~,ic] = unique(wsps,'stable');
a_counts = accumarray(ic,1);
[wsps_unique,sorted_idx_unique] = sort(wsps_unique);
[wsps_sorted,sorted_idx] = sort(wsps);
n_wsps = length(wsps_unique);
a_counts = a_counts(sorted_idx_unique);
n_seeds = max(a_counts);

seeds = NaN(n_seeds,n_wsps);
for j = 1:n_wsps
    for i = 1:a_counts(j)
        a = raw.File(ic==j);
        seeds(i,j) = str2double(regexp(a{i},'\d+(?=.dat)','match'));
    end
end

seeds(:,sorted_idx_unique);
[~,sorted_seed_idx] = sort(seeds,1);

B = reshape(sorted_idx,n_seeds,[]);
for i = 1:n_wsps
    B(:,i) = B(sorted_seed_idx(:,i),i);
end

files = cell(n_seeds,n_wsps);
data  = NaN(n_seeds,n_wsps,n_channels);
for j = 1:n_wsps
    for i = 1:n_seeds
        files{i,j} = raw.File{B(i,j),1};
        for k = 1:n_channels
            data(i,j,k) = raw{B(i,j),k+1};
        end
    end
end
end