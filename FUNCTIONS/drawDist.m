function drawDist(fig, edges, dat, color)
    figure(fig); hold on;
    histogram(dat, edges, 'Normalization', 'probability', 'FaceColor',color, 'FaceAlpha', 0.3, 'DisplayStyle','bar', 'EdgeColor','none');
    histogram(dat, edges, 'Normalization', 'probability', 'DisplayStyle','stairs', 'EdgeColor', color,'LineWidth',3, 'EdgeAlpha', 0.8);
end