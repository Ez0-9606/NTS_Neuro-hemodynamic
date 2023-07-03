function drawCumErrFrac(fig, edges, error, color)
    figure(fig); hold on;
    bin = 0.5*(edges(2:end) + edges(1:end-1));
    

    if size(error,2) > 1
        for i = 1:size(error,2)
            drawCumErrFracLight(fig, edges, error(:,i), color);
        end
    end

    [cerrTmp, ~] = histcounts(error(:), edges);
    cerr = cumsum(cerrTmp/sum(cerrTmp));
    plot(bin, cerr, 'LineWidth', 3, 'Color', color);
end