project(ngstents)

cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. $ENV{NETGENDIR}/../Resources/CMake /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen $ENV{NGSolve_DIR}
)

function(add_ngsolve_python_module2 target)
    if(NOT CMAKE_BUILD_TYPE)
      message(STATUS "Setting build type to NGSolve build type: ${NGSOLVE_CMAKE_BUILD_TYPE}")
      set(CMAKE_BUILD_TYPE ${NGSOLVE_CMAKE_BUILD_TYPE} PARENT_SCOPE)
    endif()
    add_library(${target} SHARED ${ARGN})

    #find_package (Python3 COMPONENTS Interpreter Development)
    #target_include_directories(${target} PRIVATE ${Python3_INCLUDE_DIRS})
    find_package(PythonInterp 3 REQUIRED)
    find_package(PythonLibs 3 REQUIRED)
    target_include_directories(${target} PRIVATE ${PYTHON_INCLUDE_DIR})

    if(NETGEN_BUILD_FOR_CONDA)
        message(STATUS "netgen build for conda")
        if(APPLE)
            target_link_options(${target} PUBLIC -undefined dynamic_lookup)
        else(APPLE)
        #target_link_libraries(${target} PUBLIC ${Python3_LIBRARIES})
        target_link_libraries(${target} PUBLIC ${PYTHON_LIBRARY})
        endif(APPLE)
    endif(NETGEN_BUILD_FOR_CONDA)

    set_target_properties(${target} PROPERTIES PREFIX "" CXX_STANDARD 17)

    #if(NETGEN_USE_GUI AND WIN32)
    target_link_libraries(${target} PUBLIC ngsolve)
    #else(NETGEN_USE_GUI AND WIN32)
    #  target_link_libraries(${target} PUBLIC solve)
    #endif(NETGEN_USE_GUI AND WIN32)

    if(WIN32)
      set_target_properties( ${target} PROPERTIES SUFFIX ".pyd" )
    else(WIN32)
      set_target_properties(${target} PROPERTIES SUFFIX ".so")
    endif(WIN32)

    set_target_properties(${target} PROPERTIES INSTALL_RPATH "${NGSOLVE_LIBRARY_DIR}")
endfunction()

find_package (Python3 COMPONENTS Interpreter Development)

add_ngsolve_python_module2(_pytents
  python_tents.cpp
  tents.cpp
  )
add_ngsolve_python_module2(_pyconslaw
  python_conslaw.cpp
  burgers.cpp
  euler.cpp
  wave.cpp
  advection.cpp
  maxwell.cpp
  symbolic.cpp
  vis3d.cpp
  )
target_link_libraries(_pyconslaw PRIVATE _pytents)
target_compile_definitions(_pytents PRIVATE NGSTENT_EXPORTS)

# check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR}/${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

option(BUILD_STUB_FILES "Generate stub files for better autocompletion support" OFF)
if(BUILD_STUB_FILES)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import pybind11_stubgen; print(pybind11_stubgen.__file__)" OUTPUT_VARIABLE stubgen_path RESULT_VARIABLE pybind11_stubgen)
  if(pybind11_stubgen AND NOT ${pybind11_stubgen} EQUAL 0)    
    message(WARNING "pybind11-stubgeN not found, if you want to create stub files
for better autocompletion support install it with pip.")
  else()
    message("-- Found pybind11-stubgen: ${stubgen_path}")
    set(FOUND_STUBGEN ON)
  endif()
   message(STATUS
      "This project is only compatible with\n\t\tpybind11-stubgen==0.5")
endif(BUILD_STUB_FILES)

message("With 'make install' the python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS _pytents DESTINATION ngstents)
install(TARGETS _pyconslaw DESTINATION ngstents/conslaw)


install(FILES
  ../py/__init__.py
  DESTINATION ngstents)

install(FILES
  ../py/utils/__init__.py
  ../py/utils/_drawtents.py
  ../py/utils/_drawtents2d.py
  DESTINATION ngstents/utils)

install(FILES
  ../py/conslaw/__init__.py
  DESTINATION ngstents/conslaw)

if(FOUND_STUBGEN)
  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pybind11_stubgen --no-setup-py  ngstents)")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/stubs/ngstents-stubs/ DESTINATION ngstents)
endif(FOUND_STUBGEN)
