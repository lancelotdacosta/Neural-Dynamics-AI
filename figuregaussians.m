%Figure 2
mu(1) = -2;
sigma(1) = 1;
X=-10:0.0001:10;
Y1 = normpdf(X,mu(1),sigma(1));
mu(2) = 2;
Y2 = normpdf(X,mu(2),sigma(1));
col= linspecer(4);
a=area(X,Y1,'FaceColor',col(1,:), 'FaceAlpha', 0.7);
hold on
b=area(X,Y2,'FaceColor',col(2,:), 'FaceAlpha', 0.7);
hold off
drawnow

figure(2)
sigma(2) = 4;
X=-10:0.0001:10;
Y1 = normpdf(X,mu(1),sigma(2));
Y2 = normpdf(X,mu(2),sigma(2));
a=area(X,Y1,'FaceColor',col(3,:), 'FaceAlpha', 0.7);
hold on
b=area(X,Y2,'FaceColor',col(4,:), 'FaceAlpha', 0.7);
hold off
drawnow

figure(3)
c=col;
x = [mu(1), mu(2),mu(1),   mu(2)];
y= [sigma(1), sigma(1),sigma(2),    sigma(2)];
scatter(x,y, 300, c, 'filled')
hold on
ylim([0 5])
xlim([-4 4])

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

axp = get(gca,'Position');

% determine startpoint and endpoint for the arrows 
xs=axp(1);
xe=axp(1)+axp(3)+0.04;
ys=axp(2);
ye=axp(2)+axp(4)+0.05;

% make the arrows
annotation('arrow', [xs xe],[ys ys]);
annotation('arrow', [(xs+xe-0.04+0.001)/2 (xs+xe-0.04+0.001)/2],[ys ye]);
set(gca,'color',[0.7 0.7 0.7])
hold off

drawnow




%Surface plot
% s = surf(peaks);
% alpha(s,'z')