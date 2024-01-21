function  [k_margin,d_margin,Y_margin,X_margin] = margin_adjustment(k,d,Na_C,Nr_C,Na,Nr)
[~,~,N] = size(k);
if(Na_C(end) ~= Na && Nr_C(end) ~= Nr)
    Y_margin = [1 Na_C Na]'*ones(1,length(Nr_C)+2);
    X_margin = ones(length(Na_C)+2,1)*[1 Nr_C Nr];
    
    k_margin = zeros(length(Na_C)+2,length(Nr_C)+2,N);
    k_margin(2:end-1,:,:) = [k(:,1,:) k k(:,end,:)];
    k_margin(1,:,:) = k_margin(2,:,:);
    k_margin(end,:,:) = k_margin(end-1,:,:);
    
    d_margin = zeros(length(Na_C)+2,length(Nr_C)+2,N);
    d_margin(2:end-1,:,:) = [d(:,1,:) d d(:,end,:)];
    d_margin(1,:,:) = d_margin(2,:,:);
    d_margin(end,:,:) = d_margin(end-1,:,:);
elseif(Na_C(end) == Na && Nr_C(end) == Nr)
    Y_margin = [1 Na_C]'*ones(1,length(Nr_C)+1);
    X_margin = ones(length(Na_C)+1,1)*[1 Nr_C];
    
    k_margin = zeros(length(Na_C)+1,length(Nr_C)+1,N);
    k_margin(2:end,:,:) = [k(:,1,:) k];
    k_margin(1,:,:) = k_margin(2,:,:);
    
    d_margin = zeros(length(Na_C)+1,length(Nr_C)+1,N);
    d_margin(2:end,:,:) = [d(:,1,:) d];
    d_margin(1,:,:) = d_margin(2,:,:);
elseif(Na_C(end) ~= Na && Nr_C(end) == Nr)
    Y_margin = [1 Na_C Na]'*ones(1,length(Nr_C)+1);
    X_margin = ones(length(Na_C)+2,1)*[1 Nr_C];
    
    k_margin = zeros(length(Na_C)+2,length(Nr_C)+1,N);
    k_margin(2:end-1,:,:) = [k(:,1,:) k];
    k_margin(1,:,:) = k_margin(2,:,:);
    k_margin(end,:,:) = k_margin(end-1,:,:);
    
    d_margin = zeros(length(Na_C)+2,length(Nr_C)+1,N);
    d_margin(2:end,:,:) = [d(:,1,:) d];
    d_margin(1,:,:) = d_margin(2,:,:);
    d_margin(end,:,:) = d_margin(end-1,:,:);
elseif(Na_C(end) == Na && Nr_C(end) ~= Nr)
    Y_margin = [1 Na_C]'*ones(1,length(Nr_C)+2);
    X_margin = ones(length(Na_C)+1,1)*[1 Nr_C];
    
    k_margin = zeros(length(Na_C)+1,length(Nr_C)+2,N);
    k_margin(2:end,:,:) = [k(:,1,:) k k(:,end,:)];
    k_margin(1,:,:) = k_margin(2,:,:);
    
    d_margin = zeros(length(Na_C)+1,length(Nr_C)+2,N);
    d_margin(2:end,:,:) = [d(:,1,:) d d(:,end,:)];
    d_margin(1,:,:) = d_margin(2,:,:);
end
end