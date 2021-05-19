function [index_onset, index_stable] = A_onset_stable(time, data_raw, spm, ...
    onset_period, stable_period, stable_threshold)

    % spm = smapling per minutes
    N = length(time);
   
    % onset time
    for i = 1:N
       if std(data_raw(i:i+spm*onset_period)) < 30 && mean(data_raw((i:i+spm*onset_period)))> 2
        index_onset = i+spm*onset_period;
        break;
       end
    end
    
    i_s = spm*stable_period;
    data = smoothdata(data_raw, 'movmedian',60);

    % stable time
    for i = index_onset:N
        b=i;
        if b+i_s > N
           index_stable = nan;
           break; 
        end
        x = time(b:b+i_s);
        y = data(b:b+i_s);
       a = polyfit(x,y,1);
       if abs(a(1)) <= stable_threshold && mean(y) > 2 && std(y) < 10 && mean(y) < 150
           index_stable = i;
           a(1)
           break;
       end
    end
    
end