q20 = importdata('q20-9.txt');
std = importdata('std-9.txt');
X = [1 : size(q20)];
figure(1);

plot(X, q20, X, std);

legend('Q_s = Q_p = 20', 'Standard', 'Location', 'Best');
xlabel('Wave number');
ylabel('Amplitude');
grid minor
set(gca,'FontSize',14)