from setuptools import setup
import os

shaders = [filename for filename in os.listdir("src/shaders")]

setup(name="tentsoglgui",
      version="0.1",
      description="Tents",
      packages=["tentsoglgui", "tentsoglgui.shader"],
      package_dir = {"tentsoglgui" : "src",
                     "tentsoglgui.shader" : "src/shaders"},
      install_requires=["ngsgui"],
      package_data = { "tentsoglgui" : [],
                       "tentsoglgui.shader" : shaders },
      include_package_data = True,
      entry_points = {"ngsgui.plugin" :
                      "tentsoglgui=tentsoglgui.tentplugin:loadPlugin"})
