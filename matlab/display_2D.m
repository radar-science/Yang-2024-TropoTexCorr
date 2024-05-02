function [] = display_2D(data,mask,Date,n1,n2,lim)
load colormapNew.mat
[~,~,N] = size(data);
figure;
n = 0;
for i = 1:n1
    for j = 1:n2
        n = n+1;
        if (n>N)
            break;
        end
        tmp = data(:,:,n);
        subplot(4,7,n);h = imagesc(wrap(tmp,lim));colormap(colormapNew);set(h,'alphadata',~isnan(mask));
        title(datestr(Date(n)))
    end
end
end