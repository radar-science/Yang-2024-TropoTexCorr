function [k_LLF,d_LLF,k_HTC] = slope_estimation(phase_ts_deramp,dem,MASK,W,r,N_ref)
w1 = 9;w2 = 13; % filtering parameters for obtaining high-frequency texture
res_step = 0.5;
Iteration = 10;
range = 40;
step = 0.0001;
% Input: phase time series; window size; overlap ratio
% Output: slope and intercept values in spatially discrete distributions

[Na,Nr,N] = size(phase_ts_deramp);                                         % number of images
w = (W-1)/2;                                                               % half of window size (square window: 2*w+1)
overlap = round(r*(2*w+1));                                                % overlap pixels

Na_C = w+1:(2*w-overlap):Na;
Nr_C = w+1:(2*w-overlap):Nr;
k_LLF = zeros(length(Na_C),length(Nr_C),N);                                % LLF: slope values in spatially discrete distributions
d_LLF = zeros(length(Na_C),length(Nr_C),N);                                % LLF: intercept values in spatially discrete distributions
k_HTC = zeros(length(Na_C),length(Nr_C),N);                                % HTC: slope values in spatially discrete distributions

for na = 1:length(Na_C)
    for nr = 1:length(Nr_C)

        % patch
        Na_c = Na_C(na);
        Nr_c = Nr_C(nr);
        u = Na_c-w; U = (u<1)*1+(u>=1)*u;
        d = Na_c+w; D = (d>Na)*Na+(d<=Na)*d;
        l = Nr_c-w; L = (l<1)*1+(l>=1)*l;
        r = Nr_c+w; R = (r>Nr)*Nr+(r<=Nr)*r;
        if (na == length(Na_C))
            D = Na;
        end
        if (nr == length(Nr_C))
            R = Nr;
        end
        tmp = nan(Na,Nr);tmp(U:D,L:R) = 1;
        mask_process = MASK;
        mask_process(isnan(tmp)) = NaN;
        
        % slove
        Results_Compare = zeros(4,N);                                      % line1-4: k_LLF; k_ILLF(iteration); k_HTC; d_LLF
        for n = 1:N
            if (n == N_ref)
                continue;
            end
            % ----------- initial value -----------%
            % mask
            mask_std = ones(Na,Nr).*mask_process;

            % iterative linear fitting
            % res_step = 0.5; % step [rad]
            % Iteration = 10; % 
            for iteration = 1:Iteration
                % linear fitting
                phase_tmp = phase_ts_deramp(:,:,n).*mask_std;
                dem_tmp = dem.*mask_std;
                coe = polyfit(dem(~isnan(mask_std)),phase_tmp(~isnan(mask_std)),1);
                cor_tmp = phase_tmp-coe(1)*dem_tmp-coe(2);
                
                % result recording
                if (iteration == 1)
                    Results_Compare(1,n) = coe(1);
                    Results_Compare(4,n) = coe(2);
                end
                
                % mask updating
                max_tmp = max(max(abs(cor_tmp)));
                max_tmp = (max_tmp>res_step)*(max_tmp-res_step)-(max_tmp<=res_step)*max_tmp;
                mask_std(abs(cor_tmp)>max_tmp) = NaN;
                if ( sum(sum(~isnan(mask_std))) < sum(sum(~isnan(mask_process)))/10 )
                    break;
                end
            end

            % ----------- texture correlation -----------%
            mask_tmp = MASK(U:D,L:R); 

            A = dem(U:D,L:R);
            A_LP = imgaussfilt(A,w1,'FilterSize',w2);A = A-A_LP;
            A_line = A(~isnan(mask_tmp)); A_line = (A_line/norm(A_line));

            % range = 40;
            % step = 0.0001;
            % left
            k_left = -range*step+coe(1);
            phase_ts_Scor_tmp = phase_ts_deramp(:,:,n)-k_left*dem;
            C = phase_ts_Scor_tmp(U:D,L:R); C(isnan(C)) = 0;
            C_LP = imgaussfilt(C,w1,'FilterSize',w2);C = C-C_LP;
            C_line = C(~isnan(mask_tmp)); C_line = (C_line/norm(C_line));
            conv_AC = sum(A_line.*C_line);
            record_left = abs(conv_AC);
            % right
            k_right = range*step+coe(1);
            phase_ts_Scor_tmp = phase_ts_deramp(:,:,n)-k_right*dem;
            C = phase_ts_Scor_tmp(U:D,L:R); C(isnan(C)) = 0;
            C_LP = imgaussfilt(C,w1,'FilterSize',w2);C = C-C_LP;
            C_line = C(~isnan(mask_tmp)); C_line = (C_line/norm(C_line));
            conv_AC = sum(A_line.*C_line);
            record_right = abs(conv_AC);            

            if(record_left>=record_right)
                k = (0:range)*step+coe(1);
                record = zeros(1,length(k));
            else
                k = (-range:0)*step+coe(1);
                record = zeros(1,length(k));
            end

            for i = 1:length(k)
                phase_ts_Scor_tmp = phase_ts_deramp(:,:,n)-k(i)*dem;
                C = phase_ts_Scor_tmp(U:D,L:R); C(isnan(C)) = 0;
                C_LP = imgaussfilt(C,w1,'FilterSize',w2);C = C-C_LP;
                C_line = C(~isnan(mask_tmp)); C_line = (C_line/norm(C_line));

                conv_AC = sum(A_line.*C_line);
                record(i) = abs(conv_AC);
            end
            
            % slope
            index = find(record == min(record));
            Results_Compare(2,n) = coe(1);
            Results_Compare(3,n) = k(index);
        end

        % result recording
        k_LLF(na,nr,:) = Results_Compare(1,:);
        d_LLF(na,nr,:) = Results_Compare(4,:);
        k_HTC(na,nr,:) = Results_Compare(3,:);
        
    end
end
end