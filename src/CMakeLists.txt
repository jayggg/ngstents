cmake_minimum_required(VERSION 3.16)

project(ngstents)

include(ngsolve_addon.cmake)

add_ngsolve_addon(_pytents
  python_tents.cpp
  tents.cpp
  )

add_ngsolve_addon(_pyconslaw
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

message("With 'make install' the python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
set(install_dir ${ADDON_INSTALL_DIR_PYTHON}/ngstents)
install(TARGETS _pytents DESTINATION ${install_dir})
install(TARGETS _pyconslaw DESTINATION ${install_dir}/conslaw)

install(FILES
  ../py/__init__.py
  DESTINATION ${install_dir})

install(FILES
  ../py/utils/__init__.py
  ../py/utils/_drawtents.py
  ../py/utils/_drawtents2d.py
  DESTINATION ${install_dir}/utils)

install(FILES
  ../py/conslaw/__init__.py
  DESTINATION ${install_dir}/conslaw)

# ngsolve_generate_stub_files(ngstents)
