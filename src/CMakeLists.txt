project(ngstents)

cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. $ENV{NETGENDIR}/../Resources/CMake /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)


add_ngsolve_python_module(tents
  python_tents.cpp
  tents.cpp
  )
add_ngsolve_python_module(conslaw
  python_conslaw.cpp
  symbolicDG.cpp
  burgers.cpp
  wave.cpp
  )
target_link_libraries(conslaw PRIVATE tents)

# check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR}/${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

message("With 'make install' the python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS tents DESTINATION ngstents)
install(TARGETS conslaw DESTINATION ngstents/conslaw)
install(FILES
  ../py/__init__.py
  ../py/tentslab.py
  DESTINATION ngstents)

install(FILES
  ../py/utils/__init__.py
  ../py/utils/mesh1d.py
  DESTINATION ngstents/utils)

install(FILES
  ../py/conslaw/__init__.py
  DESTINATION ngstents/conslaw)
