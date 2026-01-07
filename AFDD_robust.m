function [Freq_cluster,Freq_cluster_center,Freq_cluster_std,...
    kmeans_idx,Freq_cluster_idx,Phi_cluster] = ...
    AFDD_robust(acc_z_mobile_segment,acc_z_fixed_segment,N_mode,nfft,fs,Phi_position)
rng(0);  % fix the randomness of kNN
window_len = 60; % unit: second
slide_len = 5; % unit: second
window_N = floor((240 - window_len) / slide_len) + 1;
Phi_window = cell(length(acc_z_mobile_segment),window_N);
Freq_window = cell(length(acc_z_mobile_segment),window_N);

% frequency identification - segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(acc_z_mobile_segment)
    for j = 1:window_N
        mobile_window = acc_z_mobile_segment{i}((1:window_len*fs)+(j-1)*slide_len*fs);
        fixed_window = acc_z_fixed_segment{i}((1:window_len*fs)+(j-1)*slide_len*fs);
        [Phi_window{i,j},Freq_window{i,j},~] = AFDD([mobile_window;fixed_window],...
            1/fs:1/fs:window_len,N_mode,...
            'PickingMethod','auto','M',nfft,'dataPlot',0,'T',1);
    end
end

% use clustering methods to identify frequencies
Freq_all = [Freq_window{:,:}]';
Freq_all_idx = 1:length(Freq_all);

% Number of clusters
kmeans_k = N_mode;
% Perform k-means clustering
[kmeans_idx, ~] = kmeans(Freq_all, kmeans_k, 'Replicates', 20);
% Create a cell array to store cluster data
Freq_cluster = cell(1, kmeans_k);
Freq_cluster_idx = cell(1, kmeans_k);
for i = 1:kmeans_k
    Freq_cluster{i} = Freq_all(kmeans_idx == i);
    Freq_cluster_idx{i} = Freq_all_idx(kmeans_idx == i);
    Freq_cluster_idx{i} = Freq_cluster_idx{i}(abs(Freq_cluster{i}-median(Freq_cluster{i}))<median(Freq_cluster{i})*0.2);  % remove some outliers from the frequency cluster
    Freq_cluster{i} = Freq_cluster{i}(abs(Freq_cluster{i}-median(Freq_cluster{i}))<median(Freq_cluster{i})*0.2);  % remove some outliers from the frequency cluster
end
% Calculate mean of each cluster
Freq_cluster_center = cellfun(@median, Freq_cluster);
% Ascend Sort clusters by mean value
[~, sortIdx] = sort(Freq_cluster_center);
Freq_cluster_idx = Freq_cluster_idx(sortIdx);
Freq_cluster = Freq_cluster(sortIdx);
% Adjust cluster indices to reflect sorted order
new_kmeans_idx = zeros(size(kmeans_idx));
for newID = 1:kmeans_k
    % Get the original cluster index that should now be assigned to newID
    originalID = sortIdx(newID);
    % Assign new cluster IDs
    new_kmeans_idx(kmeans_idx == originalID) = newID;
end
% Update cluster indices and centers
kmeans_idx = new_kmeans_idx;
Freq_cluster_center = cellfun(@median, Freq_cluster);
Freq_cluster_std = cellfun(@std, Freq_cluster);


% mode shape identification - segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi_cluster = cell(kmeans_k,1);
for kk = 1:kmeans_k
    % locate identified mode shape points based on their clustering frequency
    Phi_cluster_idx = Freq_all_idx(Freq_cluster_idx{kk}); % select Phi points of the kk-th frequency
    for i = 1:length(Phi_cluster_idx)
        n1 = ceil(Phi_cluster_idx(i)/kmeans_k); % calculate the position of phi so that we can assign phi to different position
        n2 = rem(Phi_cluster_idx(i),kmeans_k);
        if n2 == 0
            n2 = kmeans_k;
        end
        n3 = ceil(n1/length(Phi_position));
        n4 = rem(n1,length(Phi_position));
        if n4 == 0
            n4 = length(Phi_position);
        end
        Phi_cluster{kk}(i,1) = n4; % position of this identified phi
        Phi_cluster{kk}(i,2:3) = Phi_window{n4,n3}(n2,:);
        Phi_cluster{kk}(i,4) = Phi_cluster{kk}(i,2)./Phi_cluster{kk}(i,3);
    end
end

end
