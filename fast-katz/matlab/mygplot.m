function [hg,hp] = mygplot(A,xy,s)

graphlw = 0.4;
bkgcolor = 0.83*[1,1,1];
graphclr = 0.55*[1,1,1];
redclr = [0.8,0,0];


clf; hold on; axis off; set(gcf,'Color',bkgcolor);
[lx,ly]=gplot(A,xy); hg=plot(lx,ly,'k-','LineWidth',graphlw);
set(hg,'Color',graphclr);
axis square

hp=scatter(xy(:,1),xy(:,2),4,s,'Filled');
colormap(hot);
set(gcf,'InvertHardCopy','off');
add_border(xy,bkgcolor);
hold off;