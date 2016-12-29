cpu = [1   1   1 ];
gpu1 = [27   1.5   19 ];
gpu2 = [50   1.8   36];
gpu4 = [94   2   62];

sum(cpu)
sum(gpu1)
sum(gpu2)
data = [cpu; gpu1; gpu2; gpu4];

bar(data);

legend({'Commp.',  'Extract & save',  'Total'}, 'Location','Best');
title('Speedup of GPUs over CPUs', 'fontsize', 14);
set(gca, 'xticklabel', {'24-core CPUs', '1 GPU', '2 GPUs', '4 GPUs'},'fontsize',12)
%axis([0,10,0,100])
ylabel('Speedup')
grid minor
set(gca, 'FontSize', 14);
%end

