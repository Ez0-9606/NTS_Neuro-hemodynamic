function prmFileCopy(srcPath, targetPath)

sid = fopen(srcPath, "r", "n");
tmp = fread(sid, '*char')';
fclose(sid);

tid = fopen(targetPath, 'w');
fwrite(tid, tmp, '*char');
fclose(tid);

end