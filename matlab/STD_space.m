unit = 100/(4*pi/lamda); % Unit conversion from [rad] to [cm]
mask = abs(velocity)>0.45;
tmp = zeros(size(velocity));tmp(150:200,160:210) = 1;tmp(320:420,295:410) = 1;tmp(481:end,455:end) = 1;
mask = mask.*tmp;
mask(isnan(mask_coh)) = 1;
% figure;imagesc(mask)

%% 时间域STD
Std_deramp = zeros(1,N);
Std_ERA5 = zeros(1,N);
Std_GLF = zeros(1,N);
Std_LLF = zeros(1,N);
Std_HTC = zeros(1,N);
for n = 1:N
    % No deformation model removal
    % tmp = phase_ts_deramp(:,:,n)*unit;
    % Std_deramp(n) = std(tmp(mask == 0),1,'all');
    % tmp = phase_ts_ERA5(:,:,n)*unit;
    % Std_ERA5(n) = std(tmp(mask == 0),1,'all');    
    % tmp = phase_ts_GLF(:,:,n)*unit;
    % Std_GLF(n) = std(tmp(mask == 0),1,'all');
    % tmp = phase_ts_LLF(:,:,n)*unit;
    % Std_LLF(n) = std(tmp(mask == 0),1,'all');
    % tmp = phase_ts_HTC(:,:,n)*unit;
    % Std_HTC(n) = std(tmp(mask == 0),1,'all');
    
    % Removal of the second-order deformation model
    tmp = phase_ts_deramp_res(:,:,n);
    Std_deramp(n) = std(tmp(mask == 0),1,'all');
    tmp = phase_ts_ERA5_res(:,:,n);
    Std_ERA5(n) = std(tmp(mask == 0),1,'all');    
    tmp = phase_ts_GLF_res(:,:,n);
    Std_GLF(n) = std(tmp(mask == 0),1,'all');
    tmp = phase_ts_LLF_res(:,:,n);
    Std_LLF(n) = std(tmp(mask == 0),1,'all');
    tmp = phase_ts_HTC_res(:,:,n);
    Std_HTC(n) = std(tmp(mask == 0),1,'all');
end

%% Display
figure;
scatter(ones(1,N),Std_deramp,'<','filled','MarkerFaceColor',[0 0 0]);hold on
scatter(1,mean(Std_deramp),55,'^','MarkerEdgeColor','b');hold on

scatter(2*ones(1,N),Std_ERA5,'<','filled','MarkerFaceColor',[0 1 0]);hold on
scatter(2,mean(Std_ERA5),55,'^','MarkerEdgeColor','b');hold on

scatter(3*ones(1,N),Std_GLF,'<','filled','MarkerFaceColor',[166/255 166/255 166/255]);hold on
scatter(3,mean(Std_GLF),55,'^','MarkerEdgeColor','b');hold on

scatter(4*ones(1,N),Std_LLF,'<','filled','MarkerFaceColor',[0 0 1]);hold on
scatter(4,mean(Std_LLF),55,'^','MarkerEdgeColor','b');hold on

scatter(5*ones(1,N),Std_HTC,'<','filled','MarkerFaceColor',[1 0 0]);hold on
scatter(5,mean(Std_HTC),55,'^','MarkerEdgeColor','b');hold on
xlim([1 5]);
boxplot([Std_deramp',Std_ERA5',Std_GLF',Std_LLF',Std_HTC'],'Widths',0.3)
ylim([0 1.5])
