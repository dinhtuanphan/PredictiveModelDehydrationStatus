clear variables;

filename2 = 'filelist.xlsx';
[~,list] = xlsread(filename2, 'Sheet1');

for name = 1: length(list)
     
close all;
clearvars -except list name

filename = char(list(name)); 
ix = strfind(filename, '.');
figurename = filename(1:ix-1)
filename_w = strcat(figurename,'_0results.xlsx');


% Slope condition
threshold = 0.1;


%% data load
raw = xlsread(filename, 'conc', 'A:E');
SR_raw = xlsread(filename, 'SR', 'A:E');

%if ~isempty(raw)
time_raw = raw(:,1);
time = 0:1/60:40;                                                     %% Set time period
%time = 0:1/60:time_raw(end);
time = time';
C_raw = raw(:,5);
%C_interp = interp1(time_raw,C_raw,time,'linear'); 
C_smooth = smoothdata(C_raw, 'movmedian',60);

time_SR_temp = SR_raw(:,1);
SR_temp = SR_raw(:,5);

k = 1;

for i =1:length(time_SR_temp)
   if ~isnan(time_SR_temp(i)) && ~isnan(SR_temp(i))
      time_SR(k) = time_SR_temp(i);
      SR(k) = SR_temp(i);
      k = k+1;
   end
end

for i = i:length(SR)
   if SR(i) >= 0
       SR(i) = SR(i);      
   end
   if SR(i) < 0
       SR(i) = 0; 
   end
end



%% onset-time && stabilization 

[i_onset,~] = A_onset_stable(time_raw,C_raw,60,0.5,5,0.1);
[~,i_stable] = A_onset_stable(time_raw,C_smooth,60,0.5,5,1);

time_onset = time(i_onset);

if ~isnan(i_stable)
    time_stable = time(i_stable);
end

if isnan(i_stable)
   time_stable = nan;  
end

[~, i_5min] = min(abs(time - 5));
[~, i_10min] = min(abs(time -10)); [~, i_15min] = min(abs(time - 15));
[~, i_30min] = min(abs(time - 30)); [~, i_40min] = min(abs(time - 40));

%% 4th order polynomial C & SR
C_poly_coeff = nan(1,5);
C_poly = nan(length(time),1);

if isnan(i_stable)
    C_poly_coeff = polyfit(time_raw(i_5min:end),C_smooth(i_5min:end),4);
    C_poly(i_5min:end) = polyval(C_poly_coeff, time(i_5min:end));
end

if ~isnan(i_stable)
    C_poly_coeff = polyfit(time_raw(i_stable:end),C_smooth(i_stable:end),4);
    C_poly(i_stable:end) = polyval(C_poly_coeff, time(i_stable:end));
end

%C_poly(i_stable:end) = polyval(C_poly_coeff, time(i_stable:end));

[~, i_SR_poly_s] = min(abs(time -time_SR(1)));
[~, i_SR_poly_e] = min(abs(time -time_SR(end)));
SR_poly = nan(length(time),1);

SR_poly_coeff = polyfit(time_SR,SR,4);
SR_poly(i_SR_poly_s:i_SR_poly_e) = polyval(SR_poly_coeff, time(i_SR_poly_s:i_SR_poly_e));

for i = 1:length(time)
    if SR_poly(i) < 0
        SR_poly(i) = nan;
    end
end

%% Linear regression of C
C_LinearR_1st = nan(length(time),1);
C_LinearR_2nd = nan(length(time),1);
C_LinearR_total = nan(length(time),1);

if time_stable < 15
    LinearR_1st = polyfit(time(i_stable:i_15min), C_poly(i_stable:i_15min),1);
    LinearR_2nd = polyfit(time(i_15min:i_40min), C_poly(i_15min:i_40min),1);
    LinearR_total = polyfit(time(i_stable:i_40min), C_poly(i_stable:i_40min),1);

    C_LinearR_1st(i_stable:i_15min) = polyval(LinearR_1st,time(i_stable:i_15min));
    C_LinearR_2nd(i_15min:i_40min) = polyval(LinearR_2nd,time(i_15min:i_40min));
    C_LinearR_total(i_stable:i_40min) = polyval(LinearR_total,time(i_stable:i_40min));

end

if time_stable >= 15
    LinearR_1st = nan(1,2);
    LinearR_2nd = polyfit(time(i_stable:i_40min), C_poly(i_stable:i_40min),1);
    LinearR_total = polyfit(time(i_stable:i_40min), C_poly(i_stable:i_40min),1);

    C_LinearR_1st = nan;
    C_LinearR_2nd(i_stable:i_40min) =  polyval(LinearR_2nd,time(i_stable:i_40min));
    C_LinearR_total(i_stable:i_40min) = polyval(LinearR_total,time(i_stable:i_40min));
    
end

if isnan(i_stable) 
    LinearR_1st = nan(1,2);
    LinearR_2nd = nan(1,2);
    LinearR_total = nan(1,2);

end

%% Linear reegression SR
%[~, i_SR_poly_s] = min(abs(time -time_SR(1)));
%[~, i_SR_poly_e] = min(abs(time -time_SR(end)));
SR_LinearR = nan(length(time),1);

SR_LinearR_coeff = polyfit(time_SR,SR,1);
SR_LinearR(i_SR_poly_s:i_SR_poly_e) = polyval(SR_LinearR_coeff, time(i_SR_poly_s:i_SR_poly_e));

for i = 1:length(time)
    if SR_LinearR(i) < 0
        SR_LinearR(i) = nan;
    end
end

%% Correlation value between C and SR
[rho, pval] = corr(C_poly,SR_poly, 'rows', 'complete');
C_SR_correl = [rho, pval];

% Maximum and Minimum value, Mean, SD

SR_max = max(SR_poly); SR_min = min(SR_poly); SR_m = mean(SR_poly);

%% Calcualtng mean value
C_mean_1st = nan;
C_mean_2nd = nan;
C_mean_total = nan;
dC = nan;
C_max = nan;
C_min = nan;
C_SD = nan;

if ~isnan(i_stable)
    C_max = max(C_poly(i_stable:i_40min)); C_min = min(C_poly(i_stable:i_40min)); 
    C_SD = std(C_poly(i_stable:i_40min));
end

if time_stable < 15
        C_mean_1st = mean(C_poly(i_stable:i_15min)); %mean value between stable to 15 minutes
        C_mean_2nd = mean(C_poly(i_15min+1:i_40min)); %mean value between 15 to 40 minutes
        C_mean_total = mean(C_poly(i_stable:i_40min));
        dC = C_mean_1st - C_mean_2nd;
end

if time_stable >= 15
        C_mean_1st = nan;
        C_mean_2nd = mean(C_poly(i_stable:i_40min)); %mean value between 15 to 40 minutes
        C_mean_total = mean(C_poly(i_stable:i_40min));
end



%% 5 minute to end

%4th order polynomial
[~, i_5min] = min(abs(time - 5));
C_poly_5min = nan(length(time),1);

C_poly_coeff_5min = polyfit(time_raw(i_5min:end),C_smooth(i_5min:end),4);
C_poly_5min(i_5min:end) = polyval(C_poly_coeff_5min, time(i_5min:end));


C_mean_1st_5min = mean(C_poly_5min(i_5min:i_15min)); %mean value between 5 to 15 minutes
C_mean_2nd_5min = mean(C_poly_5min(i_15min+1:i_40min)); %mean value between 15 to 40 minutes
C_mean_total_5min = mean(C_poly_5min(i_5min:i_40min));
dC_5min = C_mean_1st_5min - C_mean_2nd_5min;

C_max_5min = max(C_poly_5min(i_5min:i_40min)); 
C_min_5min = min(C_poly_5min(i_5min:i_40min)); 
C_SD_5min = std(C_poly_5min(i_5min:i_40min));

%Linear regression
C_LinearR_1st_5min = nan(length(time),1);
C_LinearR_2nd_5min = nan(length(time),1);
C_LinearR_total_5min = nan(length(time),1);

LinearR_1st_5min = polyfit(time(i_5min:i_15min), C_poly_5min(i_5min:i_15min),1);
LinearR_2nd_5min = polyfit(time(i_15min:i_40min), C_poly_5min(i_15min:i_40min),1);
LinearR_total_5min = polyfit(time(i_5min:i_40min), C_poly_5min(i_5min:i_40min),1);

C_LinearR_1st_5min(i_5min:i_15min) = polyval(LinearR_1st_5min,time(i_5min:i_15min));
C_LinearR_2nd_5min(i_15min:i_40min) = polyval(LinearR_2nd_5min,time(i_15min:i_40min));
C_LinearR_total_5min(i_5min:i_40min) = polyval(LinearR_total,time(i_5min:i_40min));

% Correlation value between C_5min and SR
[rho, pval] = corr(C_poly_5min,SR_poly, 'rows', 'complete');
C_SR_correl_5min = [rho, pval];

[rho, pval] = corr(C_LinearR_total_5min,SR_LinearR, 'rows', 'complete');
C_SR_correl_5min_Linear = [rho, pval];


%% Calculating Plateau

z = time_stable;
count = 1;
time_window = 5;                                               % Set the plateau window size!!

p_is = zeros(10);
p_ie = zeros(10);

%if isnan(z)
    p_count = 0;
%end

while z <= time(end) - time_window
    if isnan(z)
       break; 
    end
    [p_s, p_e, p_true] = A_plateau(time, C_poly, z, time(end), time_window, threshold, 1);
    time(p_e)  %%%% display on command window during the simulation
    
    if p_true == 1
       p_is(count) = p_s;
       p_ie(count) = p_e;
       count = count +1;
       p_count = count -1 ;
       z =  time(p_e);
    end
    
    if p_true ==0
        z = time(p_s+1);
    end      
     
end


pleatu = zeros(p_count, 2);
p_mean = zeros(p_count, 1);

plateau = zeros(p_count, 2);
for i = 1:p_count
   plateau(i,1) = time(p_is(i));
   plateau(i,2) = time(p_ie(i));
   p_mean(i) = mean(C_poly(p_is(i):p_ie(i)));
end


%% Draw figure
fig = figure('Position', [20,20,1000,700]);

subplot(2,2,1)
plot(time_raw, C_raw, 'color', [0,0,0]+0.5); hold on;
plot(time_raw, C_smooth);
xlim([0, 40]); ylim([0, 150]); hold on;

plot(time, C_poly, 'b');
hold on;

for i = 1:p_count
    plot(time(p_is(i)), p_mean(i),'r*');
    plot(time(i_onset), C_smooth(i_onset),'bo');
    hold on;
    plot(time(p_ie(i)), p_mean(i), 'r*');
    text(time(p_is(i)), p_mean(i)+5, sprintf("%.1f mM", p_mean(i)));
    line([time(p_is(i)) time(p_ie(i))], [p_mean(i) p_mean(i)], 'Color', 'red');
end

text (5, 80, sprintf("onset=%.1f ,ts = %.1f\n, C 1st=%.1f (%.1f),C 2nd=%.1f (%.1f),\n C total=%1.f (%.1f) mM, C SD = %1.f mM",...
    time_onset,time_stable, C_mean_1st, C_mean_1st_5min, C_mean_2nd, C_mean_2nd_5min, C_mean_total, C_mean_total_5min,C_SD));
figurename2 = strcat(figurename,' C');
title(figurename2);
xlabel('time (min)'); ylabel('C (mM)');
hold off;

% subplot (2,2,2)- 
%plot(time_raw, C_raw, 'color', [0,0,0]+0.5); hold on;
subplot (2,2,2)
plot(time_raw, C_smooth, 'color', [0,0,0]+0.5); hold on;
plot(time, C_poly_5min, 'b'); hold on;
xlim([0, 40]); ylim([0, 150]); hold off;


% subplot (2,2,2) - linear regression plot
%{ 
subplot(2,2,2)
plot(time_raw, C_raw, 'color', [0,0,0]+0.5); hold on;
plot(time_raw, C_smooth);
xlim([0, 40]); ylim([0, 150]); hold on;

plot(time, C_poly, 'b'); hold on;

if ~isnan(i_stable)&& time(i_stable)<=15
    plot(time, C_LinearR_1st_5min, 'b','LineWidth',2); hold on; 
    plot(time, C_LinearR_2nd_5min, 'b','LineWidth',2); hold on;
    plot(time, C_LinearR_total_5min, 'r'); hold on;
end

if time_stable > 15
   plot(time, C_LinearR_total_5min); hold on; 
end

text (5, 80, sprintf(...
    "Slope1=%.1f ,inter1 = %.1f\n,Slope2=%.1f ,inter2 = %.1f\n,Slope total=%.1f ,inter total = %.1f\n",...
    LinearR_1st_5min(1), LinearR_1st_5min(2), LinearR_2nd_5min(1), LinearR_2nd_5min(2), LinearR_total_5min(1), LinearR_total_5min(2)));
figurename4 = strcat(figurename,' C LR');
title(figurename4);
xlabel('time (min)'); ylabel('C (mM)');
hold off;
%}

subplot(2,2,3)
plot(SR_poly, C_poly_5min);
title_C_vs_SR = strcat(figurename,'(C vs SR)');
title({title_C_vs_SR});
xlabel('SR'); ylabel('C (mM)');

subplot(2,2,4)
plot(time_SR,SR, 'r*'); hold on
plot(time,SR_poly, 'b');
plot(time,SR_LinearR, 'r');
xlabel('time (min)'); ylabel('SR');
title(strcat(figurename,' SR'));
hold off;

figurename1 = strcat(figurename,'_0total_5min.jpeg');
saveas(fig,figurename1);

%{
fig0 = figure();
plot(time_raw, C_raw, 'color', [0,0,0]+0.5); hold on;
plot(time_raw, C_smooth);
xlim([0, 40]); ylim([0, 150]); hold on;

plot(time, C_poly, 'b');
hold on;

for i = 1:p_count
    plot(time(p_is(i)), p_mean(i),'r*');
    plot(time(i_onset), C_smooth(i_onset),'bo');
    hold on;
    plot(time(p_ie(i)), p_mean(i), 'r*');
    text(time(p_is(i)), p_mean(i)+5, sprintf("%.1f mM", p_mean(i)));
    line([time(p_is(i)) time(p_ie(i))], [p_mean(i) p_mean(i)], 'Color', 'red');
end

text (5, 80, sprintf("onset=%.1f ,ts = %.1f\n, C 1st=%.1f,C 2nd=%.1f,C total=%1.f mM,\n C SD = %1.f mM",...
    time_onset,time_stable, C_mean_1st, C_mean_2nd, C_mean_total, C_SD));
figurename2 = strcat(figurename,' C');
title(figurename2);
xlabel('time (min)'); ylabel('C (mM)');
hold off;

fig1 = figure();
plot(SR_poly, C_poly);
title_C_vs_SR = strcat(figurename,'(C vs SR)');
title({title_C_vs_SR;rho;pval});
xlabel('SR'); ylabel('C (mM)');

fig2 = figure();
plot(time_SR,SR, 'r*'); hold on
plot(time,SR_poly, 'b');
plot(time,SR_LinearR, 'r');
xlabel('time (min)'); ylabel('SR');
title(strcat(figurename,' SR'));
hold off;

fig3 = figure();
plot(time_raw, C_raw, 'color', [0,0,0]+0.5); hold on;
plot(time_raw, C_smooth);
xlim([0, 40]); ylim([0, 150]); hold on;

plot(time, C_poly, 'b'); hold on;

if ~isnan(i_stable)&& time(i_stable)<=15
    plot(time, C_LinearR_1st_5min, 'b','LineWidth',2); hold on; 
    plot(time, C_LinearR_2nd_5min, 'b','LineWidth',2); hold on;
    plot(time, C_LinearR_total_5min, 'r'); hold on;
end

if time_stable > 15
   plot(time, C_LinearR_total_5min); hold on; 
end

text (5, 80, sprintf(...
    "Slope1=%.1f ,inter1 = %.1f\n,Slope2=%.1f ,inter2 = %.1f\n,Slope total=%.1f ,inter total = %.1f\n",...
    LinearR_1st_5min(1), LinearR_1st_5min(2), LinearR_2nd_5min(1), LinearR_2nd_5min(2), LinearR_total_5min(1), LinearR_total_5min(2)));
figurename4 = strcat(figurename,' C LR');
title(figurename4);
xlabel('time (min)'); ylabel('C (mM)');
hold off;

figurename1 = strcat(figurename,'_1C.jpeg');
saveas(fig0,figurename1);

figurename2 = strcat(figurename,'_4C vs SR.jpeg');
saveas(fig1,figurename2);

figurename3 = strcat(figurename,'_3SR.jpeg');
saveas(fig2,figurename3);

figurename4 = strcat(figurename,'_2LinearR.jpeg');
saveas(fig3,figurename4);
%}


%% Save Excel file
header_data = {'time','C_raw', 'C_smooth', 'C_poly', 'C_LR_1st', 'C_LR_2nd', ...
    'C_LR_total', 'time_SR', 'SR_raw', 'SR_poly', 'SR_LinearR', 'C_poly_5min'};

xlswrite(filename_w, header_data, 'data', 'A1');
xlswrite(filename_w, time, 'data', 'A2');
xlswrite(filename_w, C_raw, 'data', 'B2');
xlswrite(filename_w, C_smooth, 'data', 'C2');
xlswrite(filename_w, C_poly, 'data', 'D2');
xlswrite(filename_w, C_LinearR_1st, 'data', 'E2');
xlswrite(filename_w, C_LinearR_2nd, 'data', 'F2');
xlswrite(filename_w, C_LinearR_total, 'data', 'G2');
xlswrite(filename_w, time_SR', 'data', 'H2');
xlswrite(filename_w, SR', 'data', 'I2');
xlswrite(filename_w, SR_poly, 'data', 'J2');
xlswrite(filename_w, SR_LinearR, 'data', 'K2');
    
xlswrite(filename_w, C_poly_5min, 'data', 'L2');
    
%{
header_p = {'subject', 'p_s', 'p_e', 'p_mean', 'C_last3min','onset_time','stable time', 'C_stable_30min', 'C_stable_40min','C_s1', 'C_s2', 'dC', 'slope10to40', 'pearson', 'pearson_pval','C_poly_coeff1','C_poly_coeff2','C_poly_coeff3','C_poly_coeff4','C_poly_coeff5','SR_max', 'SR_min', 'SRm', 'SR_poly_coeff1','SR_poly_coeff2'};
xlswrite(filename_w, header_p, 'summary', 'A1');
xlswrite(filename_w, {figurename}, 'summary', 'A2');

if ~isempty(plateau)
    xlswrite(filename_w, plateau, 'summary', 'B2');
    xlswrite(filename_w, p_mean, 'summary', 'D2');
end
%xlswrite(filename_w, C_last_3min, 'summary', 'E2');
xlswrite(filename_w, time_onset, 'summary', 'F2');
xlswrite(filename_w, time(i_stable), 'summary', 'G2');
xlswrite(filename_w, C_stable_30min, 'summary', 'H2');
xlswrite(filename_w, C_stable_40min, 'summary', 'I2');
xlswrite(filename_w, C_s1, 'summary', 'J2');
xlswrite(filename_w, C_s2, 'summary', 'K2');
xlswrite(filename_w, dC, 'summary', 'L2');
xlswrite(filename_w, slope, 'summary', 'M2');
xlswrite(filename_w, rho, 'summary', 'N2');
xlswrite(filename_w, pval, 'summary', 'O2');
xlswrite(filename_w, C_poly_coeff(1), 'summary', 'P2');
xlswrite(filename_w, C_poly_coeff(2), 'summary', 'Q2');
xlswrite(filename_w, C_poly_coeff(3), 'summary', 'R2');
xlswrite(filename_w, C_poly_coeff(4), 'summary', 'S2');
xlswrite(filename_w, C_poly_coeff(5), 'summary', 'T2');
xlswrite(filename_w, SRmax, 'summary', 'U2');
xlswrite(filename_w, SRmin, 'summary', 'V2');
xlswrite(filename_w, SRm, 'summary', 'W2');
xlswrite(filename_w, SR_poly_coeff(1), 'summary', 'X2');
xlswrite(filename_w, SR_poly_coeff(2), 'summary', 'Y2');
%}

%{
C_mean_1st, 
C_mean_2nd, 
C_mean_total,
dC,

C_mean_1st_5min, 
C_mean_2nd_5min, 
C_mean_total_5min, 
dC_5min, 

C_max_5min, 
C_min_5min, 
C_SD_5min, 

C_poly_coeff,
C_poly_coeff_5min,

C_SR_correl_5min, 
C_SR_correl_5min_Linear,

LinearR_1st_5min,
LinearR_2nd_5min, 
LinearR_total_5min,



%liniear regression coeff..
%Linear regression Figure..

%}



header_p = {'subject', 'p_s', 'p_e', 'p_mean', ...
    't_onset', 'ts',...
    'C_max', 'C_min', 'C_SD', 'C_mean_1st', 'C_mean_2nd', 'C_mean_total', 'dC'...
    'C_max_5min', 'C_min_5min', 'C_SD_5min', 'C_mean_1st_5min', 'C_mean_2nd_5min', 'C_mean_total_5min', 'dC_5min'...
    'C_poly_4','C_poly_3','C_poly_2','C_poly_1','C_poly_0',... 
    'C_Slope_1st', 'C_intercept_1st', 'C_Slope_2nd', 'C_intercept_2nd', 'C_Slope_total', 'C_intercept_total',...
    'C_poly_4_5min','C_poly_3_5min','C_poly_2_5min','C_poly_1_5min','C_poly_0_5min',... 
    'C_Slope_1st_5min', 'C_intercept_1st_5min', 'C_Slope_2nd_5min', 'C_intercept_2nd_5min', 'C_Slope_total_5min', 'C_intercept_total_5min',...
    'SR_mean', 'SR_max', 'SR_min', 'SR_poly_4','SR_poly_3','SR_poly_2','SR_poly_1','SR_poly_0', 'SR_slope', 'SR_intercept',...
    'SR_C_Correl_r', 'SR_C_Correl_p', 'SR_C_Correl_r_5min', 'SR_C_Correl_p_5min'};

data_array = [time(i_onset),time_stable,...
    C_max,C_min,C_SD,C_mean_1st, C_mean_2nd, C_mean_total, dC, ...
    C_max_5min, C_min_5min, C_SD_5min, C_mean_1st_5min, C_mean_2nd_5min, C_mean_total_5min, dC_5min,  ...
    C_poly_coeff, ...
    LinearR_1st, LinearR_2nd, LinearR_total, ...
    C_poly_coeff_5min, ...
    LinearR_1st_5min, LinearR_2nd_5min, LinearR_total_5min, ...
    SR_m, SR_max, SR_min,  SR_poly_coeff, SR_LinearR_coeff,...
    C_SR_correl, C_SR_correl_5min];

xlswrite(filename_w, header_p, 'summary', 'A1'); 
xlswrite(filename_w, {figurename}, 'summary', 'A2');

if ~isempty(plateau)
    xlswrite(filename_w, plateau, 'summary', 'B2'); 
    xlswrite(filename_w, p_mean, 'summary', 'D2');
end

xlswrite(filename_w, data_array, 'summary', 'E2');



%{
header_p = {'subject', 'p_s', 'p_e', 'p_mean', 't_onset', 'ts','C_max', 'C_min', 'C_SD', ...
    'C_Slope_1st', 'C_intercept_1st', 'C_Slope_2nd', 'C_intercept_2nd', 'C_Slope_total', ...
    'C_intercept_total', 'SR_mean', 'SR_max', 'SR_min', 'SR_slope', 'SR_intercept',...
    'SR_C_Correl_r', 'SR_C_Correl_p'};

xlswrite(filename_w, header_p, 'summary', 'A1'); 
xlswrite(filename_w, {figurename}, 'summary', 'A2');
if ~isempty(plateau)
    xlswrite(filename_w, plateau, 'summary', 'B2'); 
    xlswrite(filename_w, p_mean, 'summary', 'D2');
end

xlswrite(filename_w, time(i_onset), 'summary', 'E2');
xlswrite(filename_w, time_stable, 'summary', 'F2');
xlswrite(filename_w, C_max, 'summary', 'G2'); 
xlswrite(filename_w, C_min, 'summary', 'H2');
xlswrite(filename_w, C_SD, 'summary', 'I2');
xlswrite(filename_w, LinearR_1st, 'summary', 'J2');
xlswrite(filename_w, LinearR_2nd, 'summary', 'L2');
xlswrite(filename_w, LinearR_total, 'summary', 'N2');
xlswrite(filename_w, SR_m, 'summary', 'P2');
xlswrite(filename_w, SR_max, 'summary', 'Q2');
xlswrite(filename_w, SR_min, 'summary', 'R2');
xlswrite(filename_w, SR_LinearR_coeff, 'summary', 'S2');
xlswrite(filename_w, C_SR_correl, 'summary', 'U2');

header_p = {'subject', 'p_s', 'p_e', 'p_mean', 't_onset', 'ts','C_max', 'C_min', 'C_SD', ...
    'C_Slope_1st', 'C_intercept_1st', 'C_Slope_2nd', 'C_intercept_2nd', 'C_Slope_total', ...
    'C_intercept_total', 'SR_mean', 'SR_max', 'SR_min', 'SR_slope', 'SR_intercept',...
    'SR_C_Correl_r', 'SR_C_Correl_p'};

xlswrite(filename_w, header_p, 'summary', 'A1'); 
xlswrite(filename_w, {figurename}, 'summary', 'A2');
if ~isempty(plateau)
    xlswrite(filename_w, plateau, 'summary', 'B2'); 
    xlswrite(filename_w, p_mean, 'summary', 'D2');
end

xlswrite(filename_w, time(i_onset), 'summary', 'E2');
xlswrite(filename_w, time_stable, 'summary', 'F2');
xlswrite(filename_w, C_max, 'summary', 'G2'); 
xlswrite(filename_w, C_min, 'summary', 'H2');
xlswrite(filename_w, C_SD, 'summary', 'I2');
xlswrite(filename_w, LinearR_1st, 'summary', 'J2');
xlswrite(filename_w, LinearR_2nd, 'summary', 'L2');
xlswrite(filename_w, LinearR_total, 'summary', 'N2');
xlswrite(filename_w, SR_m, 'summary', 'P2');
xlswrite(filename_w, SR_max, 'summary', 'Q2');
xlswrite(filename_w, SR_min, 'summary', 'R2');
xlswrite(filename_w, SR_LinearR_coeff, 'summary', 'S2');
xlswrite(filename_w, C_SR_correl, 'summary', 'U2');
    %}
end

%{
%% Postdose

clearvars -except list name filename figurename filename_w

%% data load
raw = xlsread(filename, 'Postdose_Clinic', 'A:D');
if ~isempty(raw)
time_raw = raw(:,1);
time = 0:1/60:40;                                                     %% Set time period
time = time';
C_raw = raw(:,4);
C_interp = interp1(time_raw,C_raw,time,'linear'); 
C_smooth = smoothdata(C_interp, 'movmedian',60);
%data_value = movmean(data_raw, 30); %data_smooth = smoothdata(data_raw, 'sgolay');

%% onset-time && stabilization 

[i_onset,~] = A_onset_stable(time_raw,C_raw,60,0.5,5,0.1);
[~,i_stable] = A_onset_stable(time,C_smooth,60,0.5,5,0.1);

time_onset = time(i_onset);
if ~isnan(i_stable)
    time_stable = time(i_stable);
end

if isnan(i_stable)
   time_stable = nan;  
end

[~, i_10min] = min(abs(time -10));
[~, i_15min] = min(abs(time - 15));
[~, i_37min] = min(abs(time - 37));
[~, i_40min] = min(abs(time - 40));

if time(end) >= 40
    slope = polyfit(time(i_10min:i_40min), C_smooth(i_10min:i_40min),1);
end

if time(end) < 40
   slope = polyfit(time(i_10min:end), C_smooth(i_10min:end),1); 
end


z = time_stable;
count = 1;
time_window = 5;                                               % Set the plateau window size!!
threshold = 0.1;
p_is = zeros(10);
p_ie = zeros(10);

if isnan(z)
    p_count = 0;
end

while z <= time(end) - time_window
    if isnan(z)
       break; 
    end
    [p_s, p_e, p_true] = A_plateau(time, C_smooth, z, time(end), time_window, threshold, 1);
    time(p_e)
    
    if p_true == 1
       p_is(count) = p_s;
       p_ie(count) = p_e;
       count = count +1;
       p_count = count -1 ;
       z =  time(p_e);
    end
    

    if p_true ==0
        z = time(p_s+1);
    end      
     
end

pleatu = zeros(p_count, 2);
p_mean = zeros(p_count, 1);

plateau = zeros(p_count, 2);
for i = 1:p_count
   plateau(i,1) = time(p_is(i));
   plateau(i,2) = time(p_ie(i));
   p_mean(i) = mean(C_smooth(p_is(i):p_ie(i)));
end

C_s1 = mean(C_smooth(i_onset:i_15min));
C_s2 = mean(C_smooth(i_15min+1:i_40min));
C_onset = mean(C_smooth(i_onset:i_onset+60));

C_last_3min = mean(C_smooth(i_37min:end));
dC = C_s1 - C_s2;

fig0 = figure();
plot(time, C_smooth);
xlim([0, 40]); ylim([0, 150]); hold on;
for i = 1:p_count
    plot(time(p_is(i)), p_mean(i),'r*');
    plot(time(i_onset), C_smooth(i_onset),'bo');
    hold on;
    plot(time(p_ie(i)), p_mean(i), 'r*');
    text(time(p_is(i)), p_mean(i)+5, sprintf("%.1f mM", p_mean(i)));
    line([time(p_is(i)) time(p_ie(i))], [p_mean(i) p_mean(i)], 'Color', 'red');
end

text (5, 80, sprintf("onset=%.1f ,Conset=%.1f\n Cs1=%1.f, Cs2=%1.f, Clast3min=%.1f mM",time_onset,C_onset,C_s1,C_s2,C_last_3min));
figurename2 = strcat(figurename,' postdose');
title(figurename2);
xlabel('time (min)'); ylabel('C (mM)');
hold off;

figurename1 = strcat(figurename,'_postdose.jpeg');
saveas(fig0,figurename1);

header_p = {'subject', 'p_s', 'p_e', 'p_mean', 'C_last3min', 'time_onset', 'C_onset', 'C_s1', 'C_s2', 'dC', 'slope10to40'};
xlswrite(filename_w, header_p, 'postdose', 'A1');
xlswrite(filename_w, {figurename}, 'postdose', 'A2');

if ~isempty(plateau)
    xlswrite(filename_w, plateau, 'postdose', 'B2');
    xlswrite(filename_w, p_mean, 'postdose', 'D2');
end

xlswrite(filename_w, C_last_3min, 'postdose', 'E2');
xlswrite(filename_w, time_onset, 'postdose', 'F2');
xlswrite(filename_w, C_onset, 'postdose', 'G2');
xlswrite(filename_w, C_s1, 'postdose', 'H2');
xlswrite(filename_w, C_s2, 'postdose', 'I2');
xlswrite(filename_w, dC, 'postdose', 'J2');
xlswrite(filename_w, slope, 'postdose', 'K2');

end
%}
%}
%end