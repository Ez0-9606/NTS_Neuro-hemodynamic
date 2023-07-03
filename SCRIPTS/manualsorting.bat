@echo off
cd %1
call activate phy2
phy kwik-gui %2
call conda deactivate
echo finish