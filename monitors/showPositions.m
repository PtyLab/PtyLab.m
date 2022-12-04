figure(99)
plot( (obj.positions(:,2) - obj.positions(1,2)) * obj.dxo * 1e6, ...
      (obj.positions(:,1) - obj.positions(1,1)) * obj.dxo * 1e6,'ok-', ...
      'MarkerFaceColor',0 * [1,1,1], 'MarkerSize',12, 'LineWidth', 2)
xlabel('\mum'), ylabel('\mum')
set(gca, 'FontSize', 25)
set(gcf, 'Color', 'w');
axis square
