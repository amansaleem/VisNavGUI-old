figure;
t_CL = es_CL.es.smthBallSpd>1; % ~isnan(es_CL.es.smthBallSpd); % es_CL.es.smthBallSpd<1;
t_OL = es_OL.es.smthBallSpd<1 & es_OL.es.smthTrajSpd<1; % ~isnan(es_OL.es.smthBallSpd); % es_OL.es.smthBallSpd<1;
t_g  =  es_g.es.smthBallSpd>1; % ~isnan(es_g.es.smthBallSpd); % es_g.es.smthBallSpd<1;

c_AA_CL = corr(es_CL.es.powA(t_CL,:),es_CL.es.powA(t_CL,:));
c_BB_CL = corr(es_CL.es.powB(t_CL,:),es_CL.es.powB(t_CL,:));
c_AB_CL = corr(es_CL.es.powA(t_CL,:),es_CL.es.powB(t_CL,:));

c_AA_g = corr(es_g.es.powA(t_g,:),es_g.es.powA(t_g,:));
c_BB_g = corr(es_g.es.powB(t_g,:),es_g.es.powB(t_g,:));
c_AB_g = corr(es_g.es.powA(t_g,:),es_g.es.powB(t_g,:));

c_AA_OL = corr(es_OL.es.powA(t_OL,:),es_OL.es.powA(t_OL,:));
c_BB_OL = corr(es_OL.es.powB(t_OL,:),es_OL.es.powB(t_OL,:));
c_AB_OL = corr(es_OL.es.powA(t_OL,:),es_OL.es.powB(t_OL,:));

f = es_CL.es.freq;
range = f>1 & f<90;

subplot(331)
imagesc(f(range),f(range),c_AA_CL(range,range)); 
xlabel('CA1')
ylabel('CA1')
subplot(334)
imagesc(f(range),f(range),c_BB_CL(range,range)); 
xlabel('V1')
ylabel('V1')
subplot(337)
imagesc(f(range),f(range),c_AB_CL(range,range)); 
xlabel('CA1')
ylabel('V1')

subplot(332)
imagesc(f(range),f(range),c_AA_OL(range,range)); 
subplot(335)
imagesc(f(range),f(range),c_BB_OL(range,range)); 
subplot(338)
imagesc(f(range),f(range),c_AB_OL(range,range)); 

subplot(333)
imagesc(f(range),f(range),c_AA_g(range,range)); 
subplot(336)
imagesc(f(range),f(range),c_BB_g(range,range)); 
subplot(339)
imagesc(f(range),f(range),c_AB_g(range,range)); 

for n = 1:9
    subplot(3,3,n)
    set(gca,'box','off','TickDir','out','fontsize',12);
    axis xy; colorbar;
    axis square;
end

subplot(331)
title('CL')
subplot(332)
title('OL')
subplot(333)
title('Gray')
