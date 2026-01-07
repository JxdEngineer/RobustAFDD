clc
clear
close all

fs = 100;
nfft = 1024*2;
window = hanning(nfft/4);
N_mode = 12; % number of modes to be identified
freq_band = [3,20];
fontsize = 10.5;

load data

% observe signals in the time and frequency domain - filter out spikes
% visualization in the time domain
figure
tiledlayout(ceil(length(acc_z_mobile_segment)/2), 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:length(acc_z_mobile_segment)
    % Hampel filter
    windowSize = 50; 
    nSigma = 3; % Outlier threshold (3 standard deviations)
    acc_z_fixed_segment{i} = hampel(acc_z_fixed_segment{i}, windowSize, nSigma);
    acc_z_mobile_segment{i} = hampel(acc_z_mobile_segment{i}, windowSize, nSigma);
    % Bandpass filter
    acc_z_fixed_segment{i} = bandpass(acc_z_fixed_segment{i},freq_band,fs,'Steepness',[0.99,0.99]);
    acc_z_mobile_segment{i} = bandpass(acc_z_mobile_segment{i},freq_band,fs,'Steepness',[0.99,0.99]);

    nexttile
    hold on
    plot(1/fs:1/fs:length(acc_z_mobile_segment{i})/fs,acc_z_mobile_segment{i},'Color','blue')
    plot(1/fs:1/fs:length(acc_z_mobile_segment{i})/fs,acc_z_fixed_segment{i},'Color','red')
    ylim([-0.25,0.25])
    xlim([0,length(acc_z_mobile_segment{i})/fs])
    grid on
    box on
    xlabel('Time (s)','FontSize',fontsize,'FontName','Times New Roman')
    ylabel('Acc. Z (m/s^2)','FontSize',fontsize,'FontName','Times New Roman')
    set(gca, 'FontSize',fontsize,'FontName','Times New Roman')
    if i == length(acc_z_mobile_segment)
        legend('mobile','fixed','Location','best','FontName','Times New Roman')
    end
    title(['Stop ' num2str(i)])
end

% visualization in the frequency domain
% figure
% tiledlayout(ceil(length(acc_z_mobile_segment)/2), 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% for i = 1:length(acc_z_mobile_segment)
%     [psd_fixed,f] = pwelch(acc_z_fixed_segment{i},window,nfft/8,nfft,fs);
%     [psd_mobile,f] = pwelch(acc_z_mobile_segment{i},window,nfft/8,nfft,fs);
%     nexttile
%     hold on
%     plot(f,db(psd_mobile),'LineWidth',1.25,'Color','blue')
%     plot(f,db(psd_fixed),'LineWidth',1.25,'Color','red')
%     grid on
%     box on
%     xlabel('Frequency (Hz)','FontSize',fontsize)
%     ylabel('PSD (dB)','FontSize',fontsize)
%     set(gca, 'FontSize',fontsize)
%     if i == length(acc_z_mobile_segment)
%         legend('mobile','fixed','Location','best','Orientation','horizontal')
%     end
%     title(['Stop ' num2str(i)],'FontSize',fontsize)
%     ylim([-150,-60])
%     xlim([0,50])
% end

Phi_position = [0:6:(6-1)*length(acc_z_fixed_segment)];

% use robust AFDD for frequency and mode shape identification
[Freq_cluster,Freq_cluster_center,Freq_cluster_std,...
    kmeans_idx,Freq_cluster_idx, Phi_cluster_freq] =...
    AFDD_robust(acc_z_mobile_segment,acc_z_fixed_segment,N_mode,nfft,fs,Phi_position);

% visualize frequency identification results %%%%%%%%%%%%%%%%%%
close all
% Create the clustering plot %%%%%%%%%
figure
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile
hold on;
clr = hsv(length(Freq_cluster_center));
for i = 1:N_mode
    h1 = scatter(Freq_cluster_idx{i}, Freq_cluster{i}, 10,...
        'filled', 'MarkerFaceColor', clr(i,:));
    h2 = plot(1:length(kmeans_idx),Freq_cluster_center(i)*ones(1,length(kmeans_idx)), ...
        'LineStyle', '--', 'LineWidth', 1.5, 'color', 'k');
end
set(gca, 'FontSize',fontsize,'FontName','Times New Roman')
ylim([freq_band(1)-1,freq_band(2)+2])
xlim([1,length(kmeans_idx)])
xlabel({'Frequency Index','(a)'},'FontSize',fontsize,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',fontsize,'FontName','Times New Roman');
box on
grid on;
hold off;
legend([h1,h2],{'cluster sample','cluster center'},'Location','nw','Orientation','horizontal')
% title('(a)','FontWeight', 'normal')
% Create the box plot %%%%%%%%%
nexttile
hold on;
% Prepare data for boxplot
box_clusterData = [];
box_labels = [];
for i = 1:N_mode
    box_clusterData = [box_clusterData; Freq_cluster{i}];
    box_labels = [box_labels; repmat(i, length(Freq_cluster{i}), 1)];
end
% Overlay identified frequencies with dots
for i = 1:N_mode
    h1 = scatter(repmat(i, length(Freq_cluster{i}), 1), Freq_cluster{i}, 10, ...
        'filled', 'MarkerFaceColor', clr(i,:),'MarkerFaceAlpha',1);
end
% Plot cluster centers as black crosses
h2 = plot(1:N_mode, Freq_cluster_center, 'kx', 'MarkerSize', 15, 'LineWidth', 1.5);
% Plot box plots
boxplot(box_clusterData, box_labels, 'Widths', 0.5);
set(gca, 'FontSize',fontsize,'FontName','Times New Roman')
ylim([freq_band(1)-1,freq_band(2)+2])
xlabel({'Order of Modes','(b)'},'FontSize',fontsize,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',fontsize,'FontName','Times New Roman');
grid on;
hold off;
legend([h1,h2],{'cluster sample','cluster center'},'Location','nw','Orientation','horizontal')
% title('(b)','FontWeight', 'normal')
% overall settings %%%%%%%%%
set(gcf,'Units','centimeters',...
    'Position',[5,5,20,10],...
    'Renderer','painters')
% sgtitle('fixed + mobile')

% visualize mode shape identification results %%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi_cluster_stop = cell(N_mode,length(Phi_position));
Phi_cluster_stop_center = zeros(N_mode,length(Phi_position));
Phi_cluster_stop_std = zeros(N_mode,length(Phi_position));
figure
tiledlayout(3,ceil(N_mode/3), 'Padding', 'compact', 'TileSpacing', 'compact');
for kk = 1:N_mode
    % calculate the median of identified mode shape points at each sensing stop
    for i = 1:length(Phi_position)
        Phi_cluster_stop{kk,i} = Phi_cluster_freq{kk}(Phi_cluster_freq{kk}(:,1)==i,4);
        % remove outliers
        outliers = isoutlier(Phi_cluster_stop{kk,i}, 'median', 'ThresholdFactor', 1);
        Phi_cluster_stop{kk,i} = Phi_cluster_stop{kk,i}(~outliers);
        Phi_cluster_stop_center(kk,i) = median(Phi_cluster_stop{kk,i});
        Phi_cluster_stop_std(kk,i) = std(Phi_cluster_stop{kk,i});
    end
    Phi_cluster_center_max = max(abs(Phi_cluster_stop_center(kk,:)));
    Phi_cluster_stop_center(kk,:) = Phi_cluster_stop_center(kk,:)/Phi_cluster_center_max;

    % plot identified mode shapes
    nexttile;
    errorbar(Phi_position,Phi_cluster_stop_center(kk,:),Phi_cluster_stop_std(kk,:),...
        '-o','LineWidth',1.25,'Color','b','MarkerSize',5)
    % scatter plot identified mode shape points 
    for i = 1:length(Phi_position)
        hold on
        scatter(Phi_position(i),Phi_cluster_stop{kk,i}/Phi_cluster_center_max, 12,...
            'Marker','x','MarkerEdgeColor','red')
    end
    set(gca, 'FontSize',fontsize,'FontName','Times New Roman')
    if kk > N_mode-ceil(N_mode/3) % only show x label for the last row
        xlabel('X (m)','FontSize',fontsize,'FontName','Times New Roman')
    end
    if rem(kk,ceil(N_mode/3)) == 1
        ylabel('Mode shape','FontSize',fontsize,'FontName','Times New Roman')
    end
    box on
    grid on
    xlim([Phi_position(1)-0.5,Phi_position(end)+0.5])
    ylim([-2,2])
    if kk == N_mode
        legend('cluster center','cluster sample')
    end
    title(sprintf('Mode %d, %.2f \\pm %.2f Hz', ...
    kk, Freq_cluster_center(kk), Freq_cluster_std(kk)), ...
    'FontSize', fontsize, 'FontName', 'Times New Roman')
end
set(gcf,'Units','centimeters',...
    'Position',[25,5,20,10],...
    'Renderer','painters')