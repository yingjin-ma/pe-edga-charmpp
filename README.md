# PE-EDGA-Charm++ 

##In this version, the refactored MPS-to-CI, SR-CAS，and (PE)-EDGA charm++ codes should work together with QCMaquis.

   1. copy FindCharm.cmake  into    [QCMaquis]/dmrg/config/
   2. copy others           into    [QCMaquis]/dmrg/applications/srcas/ or [QCMaquis]/dmrg/applications/tools/
     + Notice that the CMakeList.txt should be updated basing on CMakeLists_in_srcas_or_tool.txt)
   3. And add the following into [QCMaquis]/dmrg/CMakeLists.txt before the "Version information"

```
  set(CHARM_ROOT "/Location/of/charm/folder" CACHE STRINGS "Location of charm folder")
  find_package(Charm)
  if(CHARM_FOUND)
    include_directories(${CHARM_INCLUDE_DIR})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L${CHARM_ROOT}/lib")
    # Link executables with the charmc wrapper
    STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>" "${CHARM_COMPILER}"
           CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")
  endif()

  ######################################################################
  # Version information
  ######################################################################
```
   + then re-compile QCmaquis

## Requisite:

  + Charm++  : https://github.com/UIUC-PPL/charm
  + QCmaquis : https://reiher.ethz.ch/software/maquis.html

## Related paper
  + to be updated soon

## Corresponding author
  + yingjin.ma@sccas.cn or yingjin_ma@163.com
