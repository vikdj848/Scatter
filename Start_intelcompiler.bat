call C:\Users\vikdj848\Desktop\Anaconda3\Scripts\activate.bat

call C:\Windows\System32\cmd.exe /E:ON /V:ON /K ""C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2018.0.124\windows\bin\ipsxe-comp-vars.bat" intel64 vs2015"

exit /b 

call f2py -c -m MCinloop6 inloop6.f90
