from setuptools import setup
import os

shaders = [filename for filename in os.listdir("src/shaders")]

setup(name="tentsngsgui",
      version="0.1",
      description="Tents",
      packages=["tentsngsgui", "tentsngsgui.shader"],
      package_dir = {"tentsngsgui" : "src",
                     "tentsngsgui.shader" : "src/shaders"},
      install_requires=["ngsgui"],
      package_data = { "tentsngsgui" : [],
                       "tentsngsgui.shader" : shaders },
      include_package_data = True,
      entry_points = {"ngsgui.plugin" :
                      "tentsngsgui=tentsngsgui.tentplugin:loadPlugin"})
