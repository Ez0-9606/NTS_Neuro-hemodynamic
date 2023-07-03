@echo off
cd %1
call activate klusta
klusta sorting.prm --overwrite
call conda deactivate
echo finish