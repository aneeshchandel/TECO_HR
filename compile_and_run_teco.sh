#!/bin/bash
executable_name="test.exe"
rm ${executable_name}
gfortran -o ${executable_name} TECO_HR_Aneesh.f90 && ./${executable_name} namelist_withHR_Aneesh.nml 
