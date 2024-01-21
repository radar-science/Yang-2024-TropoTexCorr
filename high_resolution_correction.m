function [phase_ts_HTC_high,flag] = high_resolution_correction(phase_ts_HTC_low,k_interp,MASK,t,w1,w2,Na_ref,Nr_ref)
phase_ts_HTC_high = phase_ts_HTC_low;
[Na,Nr,N] = size(phase_ts_HTC_high);
flag = zeros(Na,Nr);                                                       % Marks all pixels that have undergone iterative high-resolution correction

I = 4;
for iteration = 1:I
    COE = zeros(Na,Nr);                                                    % Matrix of fitting coefficients
    flag_tmp = zeros(Na,Nr);                                               % Marks pixels that have undergone high-resolution correction in this iteration
    
    % slope refinement
    for i = 1:Na
        for j = 1:Nr
            if(isnan(MASK(i,j)))
                continue;
            end
            
            % time function fitting
            ts = reshape(phase_ts_HTC_high(i,j,:),N,1);
            k = reshape(k_interp(i,j,:),N,1);
            
            A = [t,ones(N,1)];
            coe = pinv(A)*ts;
            ts_tmp1 = ts-coe(1)*t;
            
            A = [k,t,ones(N,1)];
            coe = pinv(A)*ts;
            ts_tmp2 = ts-coe(1)*k-coe(2)*t;
            
            % index for iteration
            velocity_1 = diff(ts_tmp1)./diff(t);
            velocity_2 = diff(ts_tmp2)./diff(t);

            velocity_diff1 = diff(velocity_1);
            velocity_diff2 = diff(velocity_2);
            
            delt_t =  0.5*(t(3:end)-t(1:end-2));
            acc_1 = velocity_diff1./delt_t;
            acc_2 = velocity_diff2./delt_t;

            index1 = std(acc_1(1:end))-std(acc_2(1:end));

            if (index1>0)
                flag(i,j) = flag(i,j)+1;
                flag_tmp(i,j) = 1;
                COE(i,j) = coe(1);
            end
        end
    end
    COE = imgaussfilt(COE,5,'FilterSize',5); % **
    phase_ts_HTC_high = phase_ts_HTC_high-k_interp.*COE;
    % phase_ts_HTC_high = phase_ts_HTC_high-phase_ts_HTC_high(Na_ref,Nr_ref,:);
    
    % intercept refinement
    for n = 1:N
        tmp = phase_ts_HTC_high(:,:,n);
        tmp1 = tmp; tmp1(flag_tmp==0)=0;
        tmp_filt = imgaussfilt(tmp1,w1,'FilterSize',w2);

        tmp_filt1 = tmp_filt;tmp_filt1(flag_tmp==0)=0;
        tmp_filt = imgaussfilt(tmp_filt1,25,'FilterSize',25);  % **

        tmp = tmp-tmp_filt;
        phase_ts_HTC_high(:,:,n) = tmp;
    end
    phase_ts_HTC_high = phase_ts_HTC_high-phase_ts_HTC_high(Na_ref,Nr_ref,:);
end
end