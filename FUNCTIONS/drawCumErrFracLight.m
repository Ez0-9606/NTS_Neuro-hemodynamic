function drawCumErrFracLight(fig, edges, error, color)
    if size(error,2) > 1
        error('drawCumErrFracLight::: error should be Nx1 vector\n');
    end
    figure(fig); hold on;
    bin = 0.5*(edges(2:end) + edges(1:end-1));
    [cerrTmp, ~] = histcounts(error, edges);
    cerr = cumsum(cerrTmp/sum(cerrTmp));
    plot(bin, cerr, 'LineWidth', 1, 'Color', [color, 0.2]);
end