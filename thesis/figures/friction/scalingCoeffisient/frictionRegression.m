function [] = frictionRegression( A, B )
%GRICTIONREGRESSION Summary of this function goes here
%   Detailed explanation goes here

p = polyfit(B,A,1);
l = linspace(min(B), max(B), 100);
p2 = polyval(p,l);

c1 = [63, 95, 127]/255. ;
c2 = [95, 163,142]/255.;
c3 = [0.7, 0.7, 0.7];

if p(2)>0 
    u = sprintf('$F_T=%0.2fN + %0.2f$', p(1), p(2));
else
    u = sprintf('$F_T=%0.2fN %0.2f$', p(1), p(2));
lr = plot(l, p2, 'LineWidth', 3 , 'LineStyle', '--', 'Color', c3, 'DisplayName', u );
hold on
scatter(B, A,  70, 'MarkerEdgeColor', c1, 'MarkerFaceColor', c1);
%scatter(rN,rF, 50, 'MarkerEdgeColor', c2, 'MarkerFaceColor', c2);
%scatter(N05, F05, 70, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
xlabel('$N$~[eV/\AA]',  'FontSize', 20, 'Interpreter','LaTex');
ylabel('$F_T$~[eV/\AA]','FontSize', 20, 'Interpreter','LaTex');

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

'Done'

end

