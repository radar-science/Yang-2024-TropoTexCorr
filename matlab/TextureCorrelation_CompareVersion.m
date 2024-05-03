% Function：Heterogeneous InSAR Tropospheric Correction Based on Local Texture Correlation (comparison version)
% Author：Yang Qingyue
% Date：2023.04.20
% Remark：
% 0. Methods: ERA5 correction (ERA5); global linear fitting (GLF); local linear fitting (LLF); High-frequency texture correlation (HTC)
% 1. Filter parameters：B = imgaussfilt(A,sigma,'FilterSize',w)
%    高斯滤波器的大小指定为正奇数或由正奇数组成的二元素向量。如果指定标量，则使用方滤波器. 默认滤波器大小为 2*ceil(2*sigma)+1.
% 2. Reference point：参考点选在紊流小的地方，因为方法不具备真正矫正紊流的能力，但是参考点的选取不影响低分辨率分层延迟估计过程
% 3. Algorithm acceleration：基于纹理相关的斜率求解可以使用快速搜索法，目前算法只使用了一次

clc
clear
load colormapNew

%% Parameter Setting
% About data
N = 28;                                                                    % Number of time series SAR images
N_ref = 2;                                                                 % Reference SAR image number
Na_ref = 282; Nr_ref = 204;                                                % Coordinates of the reference point: Na_ref-Azimuth; Nr_ref: Range
lamda = 0.3284035;                                                         % Radar wavelength

% About method
W = 141;                                                                   % window size (square window, must be odd number)
r = 0.4;                                                                   % overlap ratio
sigma_slope = 7; w_slope = 7;                                              % slope filtering parameters used in function imgaussfilt
sigma_intercept = 251; w_intercept = 251;                                  % intercept filtering parameters used in function imgaussfilt

%% Load Data
% acquisition time
Date = h5read('../data/timeseries_ramp_demErr.h5','/date');
Date = str2num(cell2mat(Date(1:N)));
Date = datetime(Date,'ConvertFrom','yyyymmdd');
t = days(Date-Date(1));                                                    % unit：day

% DEM
dem0 = imread('../data/height.tif');
dem = demInterp2(double(dem0));
dem(isnan(dem)) = dem0(isnan(dem));                                        % unit：meter
% figure;imagesc(dem);colormap('jet')

% rate form mintpy
rate = h5read('../data/velocity.h5','/velocity')';                         % unit: m/y
rate = double(rate*100);                                                   % unit: cm/y
% figure;imagesc(rate,[-2 2]);colormap('jet')

% deformation mask (optional, controled by maskdef)
maskdef = 0;
if (maskdef == 1)
    mask_def = rate;
    mask_def(abs(rate)>1.5) = NaN;                                         % deformation threshold(can be adjusted): 1.5m/y
    mask_def(~isnan(mask_def)) = 1;
else
    mask_def = ones(size(dem));
end

% coherence mask
mask_coh = h5read('../data/temporalCoherence.h5','/temporalCoherence')';
mask_coh = double(mask_coh);
mask_coh(mask_coh<0.87) = NaN;                                             % coherence threshold(can be adjusted): 0.87
mask_coh(~isnan(mask_coh)) = 1;
% figure;imagesc(mask_coh)

% MASK
MASK = mask_coh;
MASK = mask_def.*MASK;
% figure;imagesc(MASK)

% phase time series after：deramping + DEM error correction + ERA5 correction
phase_ts_ERA5 = h5read('../data/timeseries_ERA5_ramp_demErr.h5','/timeseries');% unit: m
phase_ts_ERA5 = double(permute(phase_ts_ERA5,[2 1 3]));
phase_ts_ERA5 = (4*pi/lamda)*phase_ts_ERA5(:,:,1:N);                         % unit: rad
phase_ts_ERA5 = phase_ts_ERA5-phase_ts_ERA5(Na_ref,Nr_ref,:);

% phase time series after：deramping + DEM error correction
phase_ts_deramp = h5read('../data/timeseries_ramp_demErr.h5','/timeseries');% unit: m
phase_ts_deramp = double(permute(phase_ts_deramp,[2 1 3]));
phase_ts_deramp = (4*pi/lamda)*phase_ts_deramp(:,:,1:N);                   % unit: rad
phase_ts_deramp = phase_ts_deramp-phase_ts_deramp(Na_ref,Nr_ref,:);
phase_ts_deramp(isnan(phase_ts_deramp)) = 0;

clear dem0

%% Local Slope Estimation Based on Texture Correlation
[Na,Nr] = size(dem);

% window segmentation
w = (W-1)/2;                                                               % half of window size (square window: 2*w+1)
overlap = round(r*(2*w+1));                                                % overlap pixels
Na_C = w+1:(2*w-overlap):Na;
Nr_C = w+1:(2*w-overlap):Nr;

% slope estimation
[k_LLF,d_LLF,k_HTC] = slope_estimation(phase_ts_deramp,dem,MASK,W,r,N_ref);% **6

%% Slope interpolation
% grid
Y = Na_C'*ones(1,length(Nr_C));
X = ones(length(Na_C),1)*Nr_C;
[Xq,Yq] = meshgrid(1:Nr,1:Na);

% K_LLF&D_LLF
[k_LLF_margin,d_LLF_margin,Y_margin,X_margin] = margin_adjustment(k_LLF,d_LLF,Na_C,Nr_C,Na,Nr);
k_LLF_interp = zeros(size(phase_ts_deramp));
d_LLF_interp = zeros(size(phase_ts_deramp));
for n = 1:N
    k_LLF_interp(:,:,n) = interp2(X_margin,Y_margin,k_LLF_margin(:,:,n),Xq,Yq,'spline');
    d_LLF_interp(:,:,n) = interp2(X_margin,Y_margin,d_LLF_margin(:,:,n),Xq,Yq,'spline');
end

% K_HTC
k_HTC_interp = zeros(size(phase_ts_deramp));
for n = 1:N
    k_HTC_filt = imgaussfilt(k_HTC(:,:,n),sigma_slope,'FilterSize',w_slope);
    k_HTC_interp(:,:,n) = interp2(X,Y,k_HTC_filt,Xq,Yq,'spline');
end

%% Intercept filtering
% GLF
k_GLF = zeros(1,N);
phase_ts_GLF = phase_ts_deramp;
for n = 1:N
    tmp = phase_ts_deramp(:,:,n);
    coe = polyfit(dem(~isnan(MASK)),tmp(~isnan(MASK)),1);
    k_GLF(n) = coe(1);
    tmp = tmp-coe(1)*dem-coe(2);
    phase_ts_GLF(:,:,n) = tmp;
end
phase_ts_GLF = phase_ts_GLF-phase_ts_GLF(Na_ref,Nr_ref,:);

% LLF
phase_ts_LLF = phase_ts_deramp;
for n = 1:N
    tmp = phase_ts_deramp(:,:,n)-k_LLF_interp(:,:,n).*dem-d_LLF_interp(:,:,n);
    phase_ts_LLF(:,:,n) = tmp;
end
phase_ts_LLF = phase_ts_LLF-phase_ts_LLF(Na_ref,Nr_ref,:);

% HTC_low
phase_ts_HTC_low = phase_ts_deramp;
intercept = zeros(size(phase_ts_HTC_low));
for n = 1:N
    tmp = phase_ts_deramp(:,:,n)-k_HTC_interp(:,:,n).*dem;
    tmp_filt = imgaussfilt(tmp,sigma_intercept,'FilterSize',w_intercept);
    tmp = tmp-tmp_filt;
    phase_ts_HTC_low(:,:,n) = tmp;
    intercept(:,:,n) = tmp_filt;
end
phase_ts_HTC_low = phase_ts_HTC_low-phase_ts_HTC_low(Na_ref,Nr_ref,:);

%% High resolution correction
% w1 = 251; w2 = 251; % **
[phase_ts_HTC_high,flag] = high_resolution_correction(phase_ts_HTC_low,k_HTC_interp,MASK,t,sigma_intercept,w_intercept,Na_ref,Nr_ref); % **2
phase_ts_HTC = phase_ts_HTC_high;

%% Velocity
tmp  = phase_ts_HTC;
velocity = zeros(Na,Nr);
for i = 1:Na
    for j = 1:Nr
        ts = reshape(tmp(i,j,:),N,1);
        coe = polyfit(t,ts,1);
        velocity(i,j) = coe(1);
    end
end
velocity = velocity*365.25*100/(4*pi/lamda);                               % unit: cm/y

figure;
h = imagesc(wrap(velocity,1.5));colormap(colormapNew);set(h,'alphadata',~isnan(mask_coh));hold on
scatter(Nr_ref,Na_ref,'k','s','filled');axis off;

%% Time series-1D
Na_c = 372; Nr_c = 372;                                                    % Coordinates of the point to be displayed

% Original
ts  = reshape(phase_ts_deramp(Na_c,Nr_c,:),1,N)/(4*pi/lamda)*100;          % unit: cm

% ERA5
ts0 = reshape(phase_ts_ERA5(Na_c,Nr_c,:),1,N)/(4*pi/lamda)*100;            % unit: cm
coe0 = polyfit(t,ts0,1);

% LLF
ts1 = reshape(phase_ts_LLF(Na_c,Nr_c,:),1,N)/(4*pi/lamda)*100;             % unit: cm
coe1 = polyfit(t,ts1,1);

% HTC
ts2 = reshape(phase_ts_HTC(Na_c,Nr_c,:),1,N)/(4*pi/lamda)*100;
coe2 = polyfit(t,ts2,1);

% slope & intercept
k = reshape(k_HTC_interp(Na_c,Nr_c,:),1,N);
d = reshape(intercept(Na_c,Nr_c,:),1,N);

figure;
subplot(2,1,1);
plot(Date,ts, ':o','MarkerSize',3,'color',[0    0    0   ],'MarkerFaceColor',[0    0    0   ],'LineWidth',1);hold on;
plot(Date,ts0,':o','MarkerSize',3,'color',[0    1    0   ],'MarkerFaceColor',[0    1    0   ],'LineWidth',1);hold on;
plot(Date,ts1,':o','MarkerSize',3,'color',[0    0    1   ],'MarkerFaceColor',[0    0    1   ],'LineWidth',1);hold on;
plot(Date,ts2,':>','MarkerSize',3,'color',[1    0    0   ],'MarkerFaceColor',[1    0    0   ],'LineWidth',1)
% plot(Date,coe0(1)*t+coe0(2),'color',[0    1    0   ])
% plot(Date,coe1(1)*t+coe1(2),'color',[0    0    1   ])
% plot(Date,coe2(1)*t+coe2(2),'color',[1    0    0   ])
% legend('Original','ERA5','Local linear fitting','The proposed method')
xlabel('Time Series');ylabel('Deformation [cm]');datetick('x','yy-mmm');ylim([-2 8])
subplot(2,1,2)
plot(Date,k);
xlabel('Time Series');ylabel('Linear cofficient');datetick('x','yy-mmm')

%% Time series-2D
[row,col] = submultiple(N);
display_2D(phase_ts_deramp,mask_coh,Date,row,col,2)                            % unit: rad
display_2D(phase_ts_ERA5,mask_coh,Date,row,col,2)                              % unit: rad
display_2D(phase_ts_GLF,mask_coh,Date,row,col,2)                               % unit: rad
display_2D(phase_ts_LLF,mask_coh,Date,row,col,2)                               % unit: rad
display_2D(phase_ts_HTC_low,mask_coh,Date,row,col,2)                           % unit: rad
display_2D(phase_ts_HTC_high,mask_coh,Date,row,col,2)                          % unit: rad

%% clear
clear coe u d l r U D L R i j na nr w w1 w2 iteration Iteration res_step initial maskdef
clear A A_line A_LP C C_line C_LP conv_AC
clear cor_tmp dem_tmp flag_tmp mask_tmp mask_std max_tmp phase_tmp phase_ts_Scor_tmp tmp tmp1 tmp_filt ts_tmp1 ts_tmp2
clear X Y Xq Yq X_line Y_line