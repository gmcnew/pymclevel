cd dist
set PATH=%PATH%;c:\program files\7-zip
7z.exe -aoa x library.zip -olibrary\
del library.zip

cd library\
7z.exe a -tzip -mx9 ..\library.zip -r
cd..
rd library /s /q

..\upx307w\upx.exe --best *.*
cd ..