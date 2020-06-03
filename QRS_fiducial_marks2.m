function [amp, ind] = QRS_fiducial_marks2(ecg_h,ecg_m,fs)

[pks,locs] = findpeaks(ecg_m,'MINPEAKDISTANCE',round(0.3*fs));

%% Step 1 - Initialize 

count = 0;           
selectedRR = 0;
meanRR = 0;        
ser_back = 0;
lengthPks = length(pks);

% Save QRS for Sig and Filtered Sig 
qrs_c = zeros(1,lengthPks);           % amplitude of R
qrs_i = zeros(1,lengthPks);           % Index
ind = zeros(1,lengthPks);             % amplitude of R
amp= zeros(1,lengthPks);              % Index

% Noise Buffers 
nois_c = zeros(1,lengthPks);
nois_i = zeros(1,lengthPks);

% Buffers for Signal and Noise
SIGL_buf = zeros(1,lengthPks);
NOISL_buf = zeros(1,lengthPks);
SIGL_buf1 = zeros(1,lengthPks);
NOISL_buf1 = zeros(1,lengthPks);
THRS_buf1 = zeros(1,lengthPks);
THRS_buf = zeros(1,lengthPks);

% Threshold signal and noise
THR_SIG = max(ecg_m(1:2*fs))*1/3;      
THR_NOISE = mean(ecg_m(1:2*fs))*1/2;   
SIG_LEV= THR_SIG;
NOISE_LEV = THR_NOISE;

% Threshold bandpass filter 
THR_SIG1 = max(ecg_h(1:2*fs))*1/3;       
THR_NOISE1 = mean(ecg_h(1:2*fs))*1/2; 
SIG_LEV1 = THR_SIG1;                     % Signal level in Band-pass filter
NOISE_LEV1 = THR_NOISE1;                 % Noise level in Band-pass filter

%% Step 2 - Threshold and decision rule

Beat_C = 0;                         % Raw Beats
Beat_C1 = 0;                        % Filtered Beats
Noise_Count = 0;                    % Noise Counter
for i = 1 : lengthPks  
    % Find the corresponding peak in the filtered signal
    if locs(i)-round(0.080*fs)>= 1 && locs(i)<= length(ecg_h)
        [y_i,x_i] = max(ecg_h(locs(i)-round(0.080*fs):locs(i)));
    else
        if i == 1
            [y_i,x_i] = max(ecg_h(1:locs(i)));
            ser_back = 1;
        elseif locs(i)>= length(ecg_h)
            [y_i,x_i] = max(ecg_h(locs(i)-round(0.080*fs):end));
        end
    end
    
    % Update the average 8RR interval
    if Beat_C >= 9 
        diffRR = diff(qrs_i(Beat_C-8:Beat_C));        % calculate RR interval
        meanRR = mean(diffRR);                        % calculate the mean of 8 previous R waves interval
        comp =qrs_i(Beat_C)-qrs_i(Beat_C-1);         
             
        if comp <= 0.92*meanRR || comp >= 1.16*meanRR   
            THR_SIG = 0.5*(THR_SIG);
            THR_SIG1 = 0.5*(THR_SIG1); 
        else
            selectedRR = meanRR;                  
        end         
    end
    
    % Identify QRS or not using average 8RR interval
    if selectedRR
        test_m = selectedRR;                        
    elseif meanRR && selectedRR == 0
        test_m = meanRR;   
    else
        test_m = 0;
    end
    
    if test_m
        if (locs(i) - qrs_i(Beat_C)) >= round(1.66*test_m)                                             
            [pks_temp,locs_temp] = max(ecg_m(qrs_i(Beat_C)+ round(0.300*fs):locs(i)-round(0.300*fs))); 
            locs_temp = qrs_i(Beat_C)+ round(0.300*fs) + locs_temp -1;      
             
            if pks_temp > THR_NOISE
                Beat_C = Beat_C + 1;
                
                qrs_c(Beat_C) = pks_temp;
                qrs_i(Beat_C) = locs_temp;      
        
                if locs_temp <= length(ecg_h)
                    [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.080*fs):locs_temp));
                else
                    [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.080*fs):end));
                end
              
                if y_i_t > THR_NOISE1 
                    Beat_C1 = Beat_C1 + 1;
                    ind(Beat_C1) = locs_temp-round(0.080*fs)+ (x_i_t - 1);
                    amp(Beat_C1) = y_i_t;                               
                    SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1;                      
                end
               
                notNoise = 1;
                SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV ;                       
            end
        else
            notNoise = 0;
        end
    end
    
    % noise or QRS peak
    if pks(i) >= THR_SIG    
        
        if Beat_C >= 3
            if (locs(i)-qrs_i(Beat_C)) <= round(0.3600*fs)
                Slope1 = mean(diff(ecg_m(locs(i)-round(0.075*fs):locs(i))));       
                Slope2 = mean(diff(ecg_m(qrs_i(Beat_C)-round(0.075*fs):qrs_i(Beat_C)))); 
                if abs(Slope1) <= abs(0.5*(Slope2))                             
                    Noise_Count = Noise_Count + 1;
                    nois_c(Noise_Count) = pks(i);
                    nois_i(Noise_Count) = locs(i);
                    count = 1;                                                 % identify T wave 
                    %noise level
                    NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
                    NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV; 
                else
                    count = 0;
                end
            end
        end
        
        % T wave detected 
        if count == 0  
            Beat_C = Beat_C + 1;
            qrs_c(Beat_C) = pks(i);
            qrs_i(Beat_C) = locs(i);
        
           if y_i >= THR_SIG1 
                Beat_C1 = Beat_C1 + 1;
                if ser_back 
                    ind(Beat_C1) = x_i;                                 
                else
                    ind(Beat_C1)= locs(i)-round(0.080*fs)+ (x_i - 1);   
                end
                amp(Beat_C1) =  y_i;                                 
                SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;                      
            end
            SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV ;                          
        end
        
    elseif (THR_NOISE <= pks(i)) && (pks(i) < THR_SIG)
         NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;                        
         NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;                            
    
    elseif pks(i) < THR_NOISE
        Noise_Count = Noise_Count + 1;
        nois_c(Noise_Count) = pks(i);
        nois_i(Noise_Count) = locs(i);    
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;                            
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;                        
    end
               
    % update threshold
    if NOISE_LEV ~= 0 || SIG_LEV ~= 0
        THR_SIG = NOISE_LEV + 0.45*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*(THR_SIG);
    end
    
    % update threshold for bandpassed signal
    if NOISE_LEV1 ~= 0 || SIG_LEV1 ~= 0
        THR_SIG1 = NOISE_LEV1 + 0.45*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*(THR_SIG1);
    end
    
    SIGL_buf(i) = SIG_LEV;
    NOISL_buf(i) = NOISE_LEV;
    THRS_buf(i) = THR_SIG;

    SIGL_buf1(i) = SIG_LEV1;
    NOISL_buf1(i) = NOISE_LEV1;
    THRS_buf1(i) = THR_SIG1;

    count = 0;                                                   
    notNoise = 0; 
    ser_back = 0;    
    
end

ind = ind(1:Beat_C1);
amp = amp(1:Beat_C1);
qrs_c = qrs_c(1:Beat_C);
qrs_i = qrs_i(1:Beat_C);

end
