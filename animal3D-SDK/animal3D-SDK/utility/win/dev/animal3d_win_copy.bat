@echo off
set srcpath=%~1
set dstpath=%~2
set srcname=%~3
set dstname=%~4
xcopy /y "%srcpath%%srcname%" "%dstpath%%dstname%"