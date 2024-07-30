close all
clear all
clc


% Load data
load RESULTS/freq.mat
load RESULTS/bayes.mat


Ea = {Ea_freq{1}, Ea_freq{2}, Ea_bayes{1}, Ea_bayes{2}, Ea_freq{3}, Ea_freq{4}, Ea_bayes{3}, Ea_bayes{4}};
lnA = {lnA_freq{1}, lnA_freq{2}, lnA_bayes{1}, lnA_bayes{2}, lnA_freq{3}, lnA_freq{4}, lnA_bayes{3}, lnA_bayes{4}};
order = [1,2,1,2,3,4,3,4];
colE_freq = [0, 0.35, 0.65];
colE_bays = [186, 218, 247]/256;
colE = {colE_freq,colE_freq, colE_bays, colE_bays, colE_freq, colE_freq, colE_bays, colE_bays };
dg = [29, 125, 80]/256;
col_dot = {'g', 'g', dg, dg, 'g', 'g', dg, dg};
edgecol = {'k', 'k', [0.8, 0.8, 0.8], [0.8, 0.8, 0.8], 'k', 'k', [0.8, 0.8, 0.8], [0.8, 0.8, 0.8]};


xE = {[-2.5, 2.5], [-2.5, 3], [-2.5, 2.5], [-2.5, 3], [-20, 60], [-30, 90], [-20, 60], [-30, 90]};

figure('Renderer', 'painters', 'Position', [200 800 1000 400])
for n = 1:8
    p = subplot(2,4,n);
    h = histogram(Ea{n});
    h.FaceColor = colE{n};
    h.EdgeColor = edgecol{n};
    h.LineWidth = 0.01;
    h.FaceAlpha = 0.95;
    p.LineWidth = 0.95;
    hold on
    scatter(mean(Ea{n}), 0, 110, col_dot{n}, 'filled')
    str = join(['Ea_', num2str(order(n))]);
    title(str, 'FontSize',17)
    if order(n)==4
        hold on
        xline(24, 'Color', 'r', 'linewidth',3)
        hold on
        xline(36, 'Color', 'r', 'linewidth',3)
    end
    if order(n)==1
        hold on
        xline(0, 'Color', 'r', 'linewidth',3)
    end
    ylim([0, 600])
    %xlim(xE{n})
end

colA_freq = [0,0.35,0.25];
colA_bays = [0.72,0.82,0.72];
colA = {colA_freq,colA_freq, colA_bays, colA_bays, colA_freq, colA_freq, colA_bays, colA_bays };

figure('Renderer', 'painters', 'Position', [200 200 1000 400])
for n = 1:8
    p = subplot(2,4,n);
    h = histogram(lnA{n});
    h.FaceColor = colA{n};
    h.EdgeColor = edgecol{n};
    h.LineWidth = 0.01;
    h.FaceAlpha = 0.95;
    p.LineWidth = 0.95;
    hold on
    scatter(mean(lnA{n}), 0, 110, col_dot{n}, 'filled')
    str = join(['ln(A_', num2str(order(n)),')']);
    title(str, 'FontSize',17)
    if order(n)==4
        hold on
        xline(log(10^13.5), 'Color', 'r', 'linewidth',4)
    end
    ylim([0, 600])
end

for r = 1:4

    mean(Ea_freq{r})

end
