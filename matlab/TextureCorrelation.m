% Function：Heterogeneous InSAR Tropospheric Correction Based on Local Texture Correlation
% Author：Yang Qingyue
% Date：2023.04.20
% Remark：
% 0. Methods: ERA5 correction (ERA5); global linear fitting (GLF); local linear fitting (LLF); High-frequency texture correlation (HTC)
% 1. Filter parameters：B = imgaussfilt(A,sigma,'FilterSize',w)
%    Size of the Gaussian filter, specified as a positive, odd integer or 2-element vector of positive, odd integers. The default filter size is 2*ceil(2*sigma)+1.
% 2. Reference point*(important)：参考点建议选在紊流小的地方，因为方法不具备真正矫正紊流的能力，但是参考点的选取不影响低分辨率分层延迟估计过程
% 3. Algorithm acceleration：To speed up the texture correlation based slope solving the fast search algorithm can be used, the algorithm is currently used only once.

clc
clear
load colormapNew

%% Parameter Setting
% About data
N = 28;                                                                    % Number of time series SAR images
N_ref = 2;                                                                 % Reference SAR image number
% Na_ref = 275; Nr_ref = 215;                                              % Coordinates of the reference point: Na_ref-Azimuth; Nr_ref: Range
Na_ref = 282; Nr_ref = 205;
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
rate = h5read('../data/velocity.h5','/velocity')';                              % unit: m/y
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
phase_ts_atm = h5read('../data/timeseries_ERA5_ramp_demErr.h5','/timeseries');% unit: m
phase_ts_atm = double(permute(phase_ts_atm,[2 1 3]));
phase_ts_atm = (4*pi/lamda)*phase_ts_atm(:,:,1:N);                         % unit: rad
phase_ts_atm = phase_ts_atm-phase_ts_atm(Na_ref,Nr_ref,:);

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
[~,~,k_HTC] = slope_estimation(phase_ts_deramp,dem,MASK,W,r,N_ref); % **6

%% Slope interpolation
% grid
Y = Na_C'*ones(1,length(Nr_C));
X = ones(length(Na_C),1)*Nr_C;
[Xq,Yq] = meshgrid(1:Nr,1:Na);

% K_HTC
k_HTC_interp = zeros(size(phase_ts_deramp));
for n = 1:N
    k_HTC_filt = imgaussfilt(k_HTC(:,:,n),sigma_slope,'FilterSize',w_slope);
    k_HTC_interp(:,:,n) = interp2(X,Y,k_HTC_filt,Xq,Yq,'spline');
end

%% Intercept filtering
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
% sigma_intercept = 251; w_intercept = 251; % **
[phase_ts_HTC_high,flag] = high_resolution_correction(phase_ts_HTC_low,k_HTC_interp,MASK,t,sigma_intercept,w_intercept,Na_ref,Nr_ref); % **2
figure;imagesc(flag);title('Number of high-resolution corrections performed')

%% Result
phase_ts_HTC = phase_ts_HTC_high;                                          % phase_ts_HTC = phase_ts_HTC_high or phase_ts_HTC_high

%% Display: Velocity
tmp  = phase_ts_HTC;                                                       % Assign tmp the data for which the deformation rate is to be calculated
% tmp = tmp-tmp(282,205,:);
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
title('Velocity [cm/y]');scatter(Nr_ref,Na_ref,'k','s','filled');axis off;

%% Display: Time series-1D
Na_c = 373; Nr_c = 373;                                                    % Coordinates of the point to be displayed

% Original
ts  = reshape(phase_ts_deramp(Na_c,Nr_c,:),1,N)/(4*pi/lamda)*100;          % unit: cm

% ERA5
ts0 = reshape(phase_ts_atm(Na_c,Nr_c,:),1,N)/(4*pi/lamda)*100;             % unit: cm
coe0 = polyfit(t,ts0,1);

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
plot(Date,ts2,':>','MarkerSize',3,'color',[1    0    0   ],'MarkerFaceColor',[1    0    0   ],'LineWidth',1)
plot(Date,coe0(1)*t+coe0(2),'color',[0    1    0   ])
plot(Date,coe2(1)*t+coe2(2),'color',[1    0    0   ])
% legend('Original','ERA5','Texture correlation')
xlabel('Time Series');ylabel('Deformation [cm]');datetick('x','yy-mmm');ylim([-2 8])
subplot(2,1,2)
plot(Date,k);
xlabel('Time Series');ylabel('Linear cofficient');datetick('x','yy-mmm')

%% Display: Time series-2D
[row,col] = submultiple(N);
display_2D(phase_ts_deramp,mask_coh,Date,row,col,2);                       % unit: rad
display_2D(phase_ts_atm,mask_coh,Date,row,col,2);                          % unit: rad
display_2D(phase_ts_HTC_low,mask_coh,Date,row,col,2);                      % unit: rad
display_2D(phase_ts_HTC_high,mask_coh,Date,row,col,2);                     % unit: rad
