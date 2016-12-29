q100 = [204 100 90 120 130 150 180 190];
q150 = [190 90 85 110 120 140 160 170];
q200 = [150 85 77 107 110 115 130 135];
q250 = [130 80 70 90 100 106 114 120];
X = [1:8];

figure(1);

plot(X, q100, '-o', X, q150, '-+', X, q200, '-*', X, q250, '-d','LineWidth',2);
xlabel('Time (s)');
ylabel('Speedup');
title('Speedup of different Q value');
legend('q=100', 'q=150', 'q=200', 'q=250', 'Location', 'Best');
grid minor
set(gca, 'FontSize', 14);

