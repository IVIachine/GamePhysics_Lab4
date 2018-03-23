@echo off
set projname=%~1
set vsver=%~2
set launchpath=%~3
set outdir=%~4
set config=%5
set buildswitch=%6
set proj="%ANIMAL3D_SDK%project\%vsver%\%projname%\%projname%.vcxproj"
::set pdb="%ANIMAL3D_SDK%%outdir%%projname%.pdb"
call "%launchpath%devenv" %proj% /%buildswitch% %config%
exit