function zach_showShift(curve1,curve2,shiftedCurve2)

figure;

a1 = axes;
a2 = axes;


plot(a1,curve1(:,1),curve1(:,2),'b','linewidth',2)
plot(a2,curve2(:,1),curve2(:,2),'r--','linewidth',2);
hold on
plot(a2,shiftedCurve2(:,1),shiftedCurve2(:,2),'r','linewidth',2);

set(a2,'color','none');
set(a2,'yaxislocation','right')
set(a2,'ycolor',[1,0,0])

set(a1,'ycolor',[0,0,1])

lim1 = get(a1,'xlim');
lim2 = get(a2,'xlim');

lim = [min([lim1(1),lim2(1)]),max([lim1(2),lim2(2)])];
set(a1,'xlim',lim)
set(a2,'xlim',lim)
set(a1,'xticklabel',[])


end