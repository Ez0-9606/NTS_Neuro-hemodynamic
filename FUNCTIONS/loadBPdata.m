function bp = loadBPdata(rawDir)

bpAddrs = findBPfile(rawDir);
bpAddrs = bpAddrs.addrs;
data = load(bpAddrs);
bp = double(data.b1(:,1));

end