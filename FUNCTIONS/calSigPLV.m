function plv = calSigPLV(firingRate, refSig)

x1 = angle(hilbert(refSig));
x2 = angle(hilbert(firingRate));
e = exp(1i*(x1-x2));
plv = abs(sum(e))/numel(x1);

end