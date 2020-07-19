from tentsngsgui import shader
from tentsngsgui.tentscene import TentScene, TentSceneSlab
from ngstents import TentSlab
import ngsgui as G


def loadPlugin(gui):
    G.shader.locations += shader.__path__
    G.GUI.sceneCreators[TentSlab] = TentScene
    G.GUI.sceneCreators[TentSlab] = TentSceneSlab
