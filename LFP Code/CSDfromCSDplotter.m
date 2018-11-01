% d0 = zeros(size(data1,1)+2,size(data1,2));
% d0(2:end-1,:) = data1;
d0 = data1(:,1:3000);
delay = 345;
% d1 = 1/3 * (diff(d0) + diff(circshift(d0,[1 0])) + diff(circshift(d0,[-1 0])));
% d2 = 1/3 * (diff(d1) + diff(circshift(d1,[1 0])) + diff(circshift(d1,[-1 0])));
% 
% % d2 = d2(3:end-2,:);
% d2c = zeros(size(d0));
% d2c(4:end-3,:) = d2(3:end-2,:);

d2c = LFP2iCSD(d0);
subplot(121)
imagesc((1:size(d0,2))./15000,1:size(d0,1), d2c);
lims = min(abs(get(gca,'Clim'))); 
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
set(gca, 'clim', [-lims lims]);
axis square;
colorbar;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
colormap(jet);

subplot(122)
plot(d2c(:,delay),(16/200):(16/200):16)
[~,pos] = min(d2c(:,delay));
axis tight
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none', 'YDir','reverse');
title(['Sink at ' num2str(16/200*pos)])
axis square; colorbar
% colorbar off
