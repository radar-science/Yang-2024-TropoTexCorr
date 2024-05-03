unit = 100/(4*pi/lamda); % Unit conversion from [rad] to [cm]

%% Original
std_deramp = zeros(Na,Nr);
phase_ts_deramp_res = zeros(Na,Nr,N);
for i = 1:Na
    for j = 1:Nr
        ts = reshape(phase_ts_deramp(i,j,:),N,1)*unit;
        coe = polyfit(t,ts,2);
        ts = ts-(coe(1)*t.^2+coe(2)*t+coe(3));
        std_deramp(i,j) = std(ts);
        phase_ts_deramp_res(i,j,:) = ts;
    end
end

%% ERA5
std_ERA5 = zeros(Na,Nr);
phase_ts_ERA5_res = zeros(Na,Nr,N);
for i = 1:Na
    for j = 1:Nr
        ts = reshape(phase_ts_ERA5(i,j,:),N,1)*unit;
        coe = polyfit(t,ts,2);
        ts = ts-(coe(1)*t.^2+coe(2)*t+coe(3));
        std_ERA5(i,j) = std(ts);
        phase_ts_ERA5_res(i,j,:) = ts;
    end
end

%% Global Linear Fitting
std_GLF = zeros(Na,Nr);
phase_ts_GLF_res = zeros(Na,Nr,N);
for i = 1:Na
    for j = 1:Nr
        ts = reshape(phase_ts_GLF(i,j,:),N,1)*unit;
        coe = polyfit(t,ts,2);
        ts = ts-(coe(1)*t.^2+coe(2)*t+coe(3));
        std_GLF(i,j) = std(ts);
        phase_ts_GLF_res(i,j,:) = ts;
    end
end

%% Local Linear fitting
std_LLF = zeros(Na,Nr);
phase_ts_LLF_res = zeros(Na,Nr,N);
for i = 1:Na
    for j = 1:Nr
        ts = reshape(phase_ts_LLF(i,j,:),N,1)*unit;
        coe = polyfit(t,ts,2);
        ts = ts-(coe(1)*t.^2+coe(2)*t+coe(3));
        std_LLF(i,j) = std(ts);
        phase_ts_LLF_res(i,j,:) = ts;
    end
end

%% Texture Correlation
std_HTC = zeros(Na,Nr);
phase_ts_HTC_res = zeros(Na,Nr,N);
for i = 1:Na
    for j = 1:Nr
        ts = reshape(phase_ts_HTC(i,j,:),N,1)*unit;
        coe = polyfit(t,ts,2);
        ts = ts-(coe(1)*t.^2+coe(2)*t+coe(3));
        std_HTC(i,j) = std(ts);
        phase_ts_HTC_res(i,j,:) = ts;
   end
end

%% Display
figure;
subplot(1,5,1);h = imagesc(std_deramp,[0 1.5]);colormap('cool');set(h,'alphadata',~isnan(mask_coh));axis off;hold on;scatter(Nr_ref,Na_ref,'k','s','filled');axis off;
subplot(1,5,2);h = imagesc(std_ERA5,[0 1.5]);colormap('cool');set(h,'alphadata',~isnan(mask_coh));axis off;hold on;scatter(Nr_ref,Na_ref,'k','s','filled');axis off;
subplot(1,5,3);h = imagesc(std_GLF,[0 1.5]);colormap('cool');set(h,'alphadata',~isnan(mask_coh));axis off;hold on;scatter(Nr_ref,Na_ref,'k','s','filled');axis off;
subplot(1,5,4);h = imagesc(std_LLF,[0 1.5]);colormap('cool');set(h,'alphadata',~isnan(mask_coh));axis off;hold on;scatter(Nr_ref,Na_ref,'k','s','filled');axis off;
subplot(1,5,5);h = imagesc(std_HTC,[0 1.5]);colormap('cool');set(h,'alphadata',~isnan(mask_coh));axis off;hold on;scatter(Nr_ref,Na_ref,'k','s','filled');axis off;

