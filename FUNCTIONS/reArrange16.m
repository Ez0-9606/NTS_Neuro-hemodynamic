function arrng = reArrange16(ntsData)
% Jiho Lee, 2022-12-12
% input:
% NTS recording data (16 = # of channel) x (# of time steps):
%
% output: re-arranged 16 channel NTS recording data
% 
arrng = zeros(size(ntsData));
Transform = getTransformMtx();
Arrangement = Transform*(1:1:16)';
for i = 1:16
    arrng(Arrangement(i),:) = ntsData(i,:);
end

end