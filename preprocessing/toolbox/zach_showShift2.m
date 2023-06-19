function f = zach_showShift2(curve1,curve2,shiftedCurve2,curve3,shiftedCurve3)

f = figure;

a1 = axes;
plot(a1,curve1(:,1),curve1(:,2),'b','linewidth',2)
hold on
xlabel('Time (s)')
ylabel('BOLD units')
set(a1,'ticklength',[0 0])

a2 = axes;
plot(a2,shiftedCurve2(:,1),shiftedCurve2(:,2),'r-','linewidth',2);
hold on
plot(a2,curve2(:,1),curve2(:,2),'r--','linewidth',0.5);
ylabel('Partial Pressure (arbitrary units)')
set(a2,'yaxislocation','right')
set(a2,'color','none');
set(a2,'ticklength',[0 0])

a3 = axes;
plot(a3,shiftedCurve3(:,1),shiftedCurve3(:,2),'m-','linewidth',2);
hold on
plot(a3,curve3(:,1),curve3(:,2),'m--','linewidth',0.5);
set(a3,'color','none');
set(a3,'ticklength',[0 0])

xlim1 = get(a1,'xlim');
xlim2 = get(a2,'xlim');
xlim3 = get(a3,'xlim');

xlim = [min([xlim1(1),xlim2(1),xlim3(1)]),max([xlim1(2),xlim2(2),xlim3(2)])];
set(a1,'xlim',xlim)
set(a2,'xlim',xlim)
set(a3,'xlim',xlim)
set(a2,'xticklabel',[])
set(a3,'xticklabel',[])

% ylim2 = get(a2,'ylim');
% ylim3 = get(a3,'ylim');
% 
% ylim = [min([ylim2(1),ylim3(1)]),max([ylim2(2),ylim3(2)])];
% set(a2,'ylim',ylim)
% set(a3,'ylim',ylim)
set(a3,'yticklabel',[])
set(a2,'yticklabel',[])




% plot(a1,shiftedCurve2(:,1),shiftedCurve2(:,2),'r-','linewidth',2);
% plot(a1,curve2(:,1),curve2(:,2),'r--','linewidth',0.5);
% 
% plot(a2,shiftedCurve3(:,1),shiftedCurve3(:,2),'m-','linewidth',2);
% plot(a2,curve3(:,1),curve3(:,2),'m--','linewidth',0.5);



% set(a1,'color','none');
% set(a1,'yaxislocation','right')
% set(a1,'ycolor',[1,0,0])
% 
% set(a1,'ycolor',[0,0,1])
% 
% lim1 = get(a1,'xlim');
% lim2 = get(a2,'xlim');
% 
% lim = [min([lim1(1),lim2(1)]),max([lim1(2),lim2(2)])];
% set(a1,'xlim',lim)
% set(a2,'xlim',lim)
% set(a1,'xticklabel',[])


end