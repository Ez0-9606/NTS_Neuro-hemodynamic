function plv = calPhPLV(phase, ref)
    plv = abs(sum(exp(1i*(phase-ref))))/numel(phase);
end