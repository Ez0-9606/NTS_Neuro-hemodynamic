function out = findBPfile(directStruct)
% Find the file name and its index of BP recording data from directory file list
out = struct();
j = 1;
    for i = 1:numel(directStruct)
        tmp = directStruct(i);
        if contains(tmp.name, ".mat")
            out(j).name = tmp.name;
            out(j).addrs = [tmp.folder, '\', tmp.name];
            out(j).idx = i;
            j = j + 1;
        end
    end

end