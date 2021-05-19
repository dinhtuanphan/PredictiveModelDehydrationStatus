clear variables;

filename2 = 'filelist.xlsx';
[~,list] = xlsread(filename2, 'Sheet1');

for name = 1: length(list)
     
close all;
clearvars -except list name

filename = char(list(name)); 
ix = strfind(filename, '.');
figurename = filename(1:ix-1)
filename_w = strcat(figurename,'_results.xlsx');


% Slope condition
threshold = 0.1;


%% data load
raw = xlsread(filename, 'conc', 'A:E');
SR_raw = xlsread(filename, 'SR', 'A:E');

%if ~isempty(raw)
time_raw = raw(:,1);
%time = 0:1/60:40;                                                     %% Set time period
time = 0:1/60:time_raw(end);
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


%% onset-time && stabilization 

[i_onset,~] = A_onset_stable(time_raw,C_raw,60,0.5,5,0.1);
[~,i_stable] = A_onset_stable(time_raw,C_smooth,60,0.5,5,1);

%% 4th order polynomial C and SR

C_poly = nan(length(time),1);
if ~isnan(i_stable)
    %C_smooth_temp = rmmissing(C_smooth);
    C_poly_coeff = polyfit(time_raw(i_stable:end),C_smooth(i_stable:end),4);
end

C_poly(i_stable:end) = polyval(C_poly_coeff, time(i_stable:end));

SR_poly_coeff = polyfit(time_SR,SR,1);

SR_poly = nan(length(time),1);
[~, i_SR_poly_s] = min(abs(time -time_SR(1)));
[~, i_SR_poly_e] = min(abs(time -time_SR(end)));

SR_poly(i_SR_poly_s:i_SR_poly_e) = polyval(SR_poly_coeff, time(i_SR_poly_s:i_SR_poly_e));

for i = 1:length(time)
    if SR_poly(i) < 0
        SR_poly(i) = nan;
    end
end

% Correlation value between C and SR
[rho, pval] = corr(C_poly,SR_poly, 'rows', 'complete');

%% slope from 10 to time(end)
time_onset = time(i_onset);

if ~isnan(i_stable)
    time_stable = time(i_stable);
end

if isnan(i_stable)
   time_stable = nan;  
end

[~, i_10min] = min(abs(time -10));
[~, i_15min] = min(abs(time - 15));
[~, i_30min] = min(abs(time - 30));
%[~, i_37min] = min(abs(time - 37));
[~, i_40min] = min(abs(time - 40));

%if time(end) >= 40
slope = polyfit(time(i_10min:i_40min), C_poly(i_10min:i_40min),1);
%end
%
%if time(end) < 40
%   slope = polyfit(time_smooth(i_10min:end), C_smooth(i_10min:end),1); 
%end

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
   p_mean(i) = mean(C_poly(p_is(i):p_ie(i)));
end

%% Calcualtng mean value

if i_stable <= i_15min
        C_s1 = mean(C_poly(i_stable:i_15min)); %mean value between stable to 15 minutes
end

if i_stable > i_15min
        C_s1 = nan;
end

C_s2 = mean(C_poly(i_15min+1:i_40min)); %mean value between 15 to 40 minutes

C_stable_30min = nan;
C_stable_40min = nan;

if ~isnan(i_stable)
    C_stable_30min = mean(C_poly(i_stable:i_30min)); %mean value between stable to 30 minutes
    C_stable_40min = mean(C_poly(i_stable:i_40min)); %mean value between stable to 40 minutes
end

%C_last_3min = mean(C_smooth(i_37min:end));
dC = C_s1 - C_s2;

%% Draw figure

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

text (5, 80, sprintf("onset=%.1f ,Cs30min=%.1f\n,Cs40min=%.1f\n  Cs1=%1.f, Cs2=%1.f mM",time_onset,C_stable_30min,C_stable_40min,C_s1,C_s2));
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
plot(time,SR_poly, 'r');
xlabel('time (min)'); ylabel('SR');
title(strcat(figurename,' SR'));
hold off;


figurename1 = strcat(figurename,'_C.jpeg');
saveas(fig0,figurename1);

figurename2 = strcat(figurename,'_C vs SR.jpeg');
saveas(fig1,figurename2);

figurename3 = strcat(figurename,'_SR.jpeg');
saveas(fig2,figurename3);


%% Save Excel file
header_data = {'time','C_raw', 'C_smooth', 'C_poly', 'SR_poly'};
xlswrite(filename_w, header_data, 'data', 'A1');
xlswrite(filename_w, time, 'data', 'A2');
xlswrite(filename_w, C_raw, 'data', 'B2');
xlswrite(filename_w, C_smooth, 'data', 'C2');
xlswrite(filename_w, C_poly, 'data', 'D2');
xlswrite(filename_w, SR_poly, 'data', 'E2');

header_p = {'subject', 'p_s', 'p_e', 'p_mean', 'C_last3min','onset_time','stable time', 'C_stable_30min', 'C_stable_40min','C_s1', 'C_s2', 'dC', 'slope10to40', 'pearson', 'pearson_pval','C_poly_coeff','SR_poly_coeff'};
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
xlswrite(filename_w, C_poly_coeff', 'summary', 'P2');
xlswrite(filename_w, SR_poly_coeff', 'summary', 'Q2');


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