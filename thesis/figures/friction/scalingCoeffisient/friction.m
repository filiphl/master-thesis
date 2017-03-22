clear all

F = [608         506         259         853         936        1400        1500         391         552         731         888        1142];
N = [1077        805         541         1341        1613       1881        2150         672         942         1208        1480       1745];


%      5        6        8       9      10      11      12       13       *14     15      16      17      18      *19     *20
F05 = [188      270     656     558     613     732     780      655      951     1064    1099    1231    1271    1483    1530]
N05 = [132      271     537     670     808     940     1079     1210     1343    1480    1615    1750    1881    2015    2151]




%rF = [1400 888  1142];
%rN = [1881 1480 1745];
p = polyfit(N,F,1);
l = linspace(min(N), max(N), 100);
p2 = polyval(p,l);

c1 = [63, 95, 127]/255. ;
c2 = [95, 163,142]/255.;
c3 = [0.7, 0.7, 0.7];

u = sprintf('$F_T=%0.2fN %0.2f$', p(1), p(2));
lr = plot(l, p2, 'LineWidth', 3 , 'LineStyle', '--', 'Color', c3, 'DisplayName', u );
hold on
scatter(N, F,  70, 'MarkerEdgeColor', c1, 'MarkerFaceColor', c1);
%scatter(rN,rF, 50, 'MarkerEdgeColor', c2, 'MarkerFaceColor', c2);
%scatter(N05, F05, 70, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
xlabel('$N$~[eV/\AA]',  'FontSize', 16, 'Interpreter','LaTex');
ylabel('$F_T$~[eV/\AA]','FontSize', 16, 'Interpreter','LaTex');
%t = sprintf('$\\mu = %0.2f$',p(1));
%title(t, 'FontSize', 18, 'interpreter','latex');

[h, hobj, plt, ~] = legend(lr);
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',3);
hl(1).XData(2) = hl(1).XData(2) - 0.05;
ht = findobj(hobj,'type','text');
set(ht,'FontSize',16);
set(ht, 'interpreter', 'latex');
set(h, 'Location', 'North west');
pos = get(ht, 'Position');
pos(1) = pos(1) - 0.02;
set(ht, 'Position', pos); 
h.Position(3) = h.Position(3)+0.16
h.Position(2) = h.Position(2)-0.02
h.Position(4) = h.Position(4)+0.02

h.EdgeColor = [0.5, 0.5, 0.5]
h.LineWidth = 1
