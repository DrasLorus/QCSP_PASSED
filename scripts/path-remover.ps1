param(
    [String]$filein,
    [String]$fileout
) 

((Get-Content $filein) -creplace '#include "\./','#include "') -creplace '#include ".*/','#include "' > $fileout