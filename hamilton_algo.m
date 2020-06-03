function [amp, ind] = hamilton_algo(ecg,fs)

%% Input data
ecg = ecg(:);    

%% Band-pass filtering(8-16Hz)

% fs=200
if fs == 200   
    
    ecg = ecg-mean(ecg);            %remove the mean of signal
    
    Wn=16*2/fs;                    
    N = 3;                         
    [a,b] = butter(N, Wn, 'low');   %band-pass filtering
    ecg_l = filtfilt(a, b, ecg);
    ecg_l = ecg_l/max(abs(ecg_l));
    
    Wn = 8*2/fs;
    N = 3;                          
    [a,b] = butter(N,Wn,'high');    % band-pass filtering
    ecg_h = filtfilt(a,b,ecg_l); 
    ecg_h = ecg_h/ max(abs(ecg_h));  
    
else  
    
    % fs ~= 200
    f1=8;                                                                      
    f2=16;                                                                     
    Wn=[f1 f2]*2/fs;                                                         
    N = 3;                                                                      
    [a,b] = butter(N,Wn);           % band-pass filtering                                         
    ecg_h = filtfilt(a,b,ecg);
    ecg_h = ecg_h/ max(abs(ecg_h));  
    
end

%% Derivative filtering

if fs ~= 200
    int_c = (5-1)/(fs*1/40);
    b = interp1(1:5,[1 2 0 -2 -1].*(1/8)*fs,1:int_c:5);
else
    b = [1 2 0 -2 -1].*(1/8)*fs;   
end

ecg_d = filtfilt(b,1,ecg_h);
ecg_d = ecg_d/max(ecg_d);

%% Rectified

ecg_s = abs(ecg_d);
 
%% Moving window average

ecg_m = conv(ecg_s ,ones(1 ,round(0.080*fs))/round(0.080*fs));% p&t 0.15*fs = 54 sample

%% Find QRS Peaks 

[amp,ind] = QRS_fiducial_marks2(ecg_h,ecg_m,fs);
 
end
