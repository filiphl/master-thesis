clear all
N = 1208;
F = [524.090, 544.118, 558.7, 653.1, 754.7, 815.8, 901.5]
u = F/N
v = [0.25 0.5 1, 2, 3, 4, 5];
UofV(u,v)



N100 = 807.382;
F100 = [612.759, 679.121, 751.159];
u100=F100/N100;
v100 = [3, 4, 5]
%UofV(u100, v100)


% p = polyfit(v,u,1);
% l = linspace(0, max(v), 100);
% p2 = polyval(p,l);
% 
% c1 = [63, 95, 127]/255. ;
% c2 = [95, 163,142]/255.;
% c3 = [0.7, 0.7, 0.7];
% if p(2)>0 
%     L = sprintf('$\\mu(v)=%0.4fv + %0.4f$', p(1), p(2))
% else
%     L = sprintf('$\\mu(v)=%0.2fv %0.2f$', p(1), p(2));
% end
% 
% 
% fig = figure;
% lr = plot(l, p2, 'LineWidth', 3 , 'LineStyle', '--', 'Color', c3, 'DisplayName', L );
% hold on
% scatter(v, u,  70, 'MarkerEdgeColor', c1, 'MarkerFaceColor', c1);
% %scatter(rN,rF, 50, 'MarkerEdgeColor', c2, 'MarkerFaceColor', c2);
% %scatter(N05, F05, 70, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
% xlabel('$v$~[\AA/ps]',  'FontSize', 18, 'Interpreter','LaTex');
% ylabel('$\frac{F_{max}}{N}$','FontSize', 22, 'Interpreter','LaTex', 'rot', 0, 'Position', [-0.7, 0.58]);
% %t = sprintf('$\\mu = %0.2f$',p(1));
% %title(t, 'FontSize', 18, 'interpreter','latex');
% 
% [h, hobj, plt, ~] = legend(lr);
% hl = findobj(hobj,'type','line');
% set(hl,'LineWidth',3);
% hl(1).XData(2) = hl(1).XData(2) - 0.05;
% ht = findobj(hobj,'type','text');
% set(ht,'FontSize',16);
% set(ht, 'interpreter', 'latex');
% set(h, 'Location', 'North west');
% pos = get(ht, 'Position');
% pos(1) = pos(1) - 0.02;
% set(ht, 'Position', pos); 
% h.Position(3) = h.Position(3)+0.16;
% h.Position(2) = h.Position(2)-0.02;
% h.Position(4) = h.Position(4)+0.02;
% 
% h.EdgeColor = [0.5, 0.5, 0.5];
% h.LineWidth = 1;
% 
% 
% set(fig,'Units','Inches');
% 
% pos = get(fig,'Position')
% 
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)*1.08, pos(4)*1.1])
