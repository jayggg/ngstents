from ngsgui.scenes import BaseMeshScene
from ngsgui.gl import Texture, getProgram
from ngsgui.thread import inmain_decorator
import ngsgui.settings as settings
from qtpy import QtCore
import numpy
import ngsolve

import numpy as np
import OpenGL.GL as GL

from enum import Enum


class DrawTents(Enum):
    ALL = 0
    TENTNR = 1
    LEVELNR = 2


class TentScene(BaseMeshScene, settings.ColormapSettings):
    @inmain_decorator(wait_for_return=True)
    def __init__(self, tentslab, ** kwargs):
        self._initial_values = {"SelectedTent": 0,
                                "ShowSingleTent": False,
                                "SelectedLevel": 0,
                                "ShowSingleLevel": False,
                                "ShrinkTents": 0.95,
                                "ScaleTents": 1.0}

        mesh = tentslab.mesh
        tentdata, tenttimes, ntents, nlevels = tentslab.GetTentData()
        self.tentdata = tentdata
        self.tenttimes = tenttimes
        self.ntents = ntents
        self.nlevels = nlevels
        self.lententdata = len(tentdata)//4

        super().__init__(mesh, **kwargs)

    def initGL(self):
        super().initGL()
        self.vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(self.vao)
        self.tex_tents = Texture(GL.GL_TEXTURE_BUFFER, GL.GL_RGBA32I)
        self.tex_times = Texture(GL.GL_TEXTURE_BUFFER, GL.GL_RGBA32F)
        self.tex_tents_color = Texture(GL.GL_TEXTURE_1D, GL.GL_RGBA)
        GL.glBindVertexArray(0)

    @inmain_decorator(True)
    def update(self):
        super().update()
        # store tent data on GPU memory
        tents = np.array(self.tentdata, dtype=np.int32)
        n = len(tents)
        self.tex_tents.store(tents, GL.GL_INT)
        tenttimes = np.zeros(n, dtype=np.float32) if self.tenttimes is None \
                    else np.array(self.tenttimes, dtype=np.float32)
        self.tex_times.store(tenttimes)
        self.tents_colors = [0, 255, 0, 255]*self.ntents
        self.tex_tents_color.store(
            self.tents_colors, GL.GL_UNSIGNED_BYTE, self.ntents)

    @inmain_decorator(True)
    def _createParameters(self):
        super()._createParameters()
        animate_tent = settings.CheckboxParameter(name="AnimateTent",
                                                  label="Animate tents")

        def animateTent(val):
            if val:
                self._timer_thread = QtCore.QThread()

                def run_animate():
                    self._animation_timer = QtCore.QTimer()
                    self._animation_timer.setInterval(20)
                    self._animation_timer.timeout.connect(
                        lambda: self.setSelectedTent(self.getSelectedTent()+1))
                    self._animation_timer.start()

                def stop_animate():
                    self._animation_timer.stop()
                self._timer_thread.started.connect(run_animate)
                self._timer_thread.finished.connect(stop_animate)
                self._timer_thread.start()
            else:
                self._timer_thread.finished.emit()
                self._timer_thread.quit()
        animate_tent.changed.connect(animateTent)

        animate_lvl = settings.CheckboxParameter(name="AnimateLevel",
                                                 label="Animate levels")

        def animateLevel(val):
            if val:
                self._timer_thread = QtCore.QThread()

                def run_animate():
                    self._animation_timer = QtCore.QTimer()
                    self._animation_timer.setInterval(100)
                    self._animation_timer.timeout.connect(
                        lambda: self.setSelectedLevel(
                            (self.getSelectedLevel()+1) % self.nlevels))
                    self._animation_timer.start()

                def stop_animate():
                    self._animation_timer.stop()
                self._timer_thread.started.connect(run_animate)
                self._timer_thread.finished.connect(stop_animate)
                self._timer_thread.start()
            else:
                self._timer_thread.finished.emit()
                self._timer_thread.quit()
        animate_lvl.changed.connect(animateLevel)

        sub_parameters = [
            [],
            [settings.ValueParameter(
                name="SelectedTent", label="Tent",
                default_value=self._initial_values["SelectedTent"],
                min_value=0, max_value=self.ntents-1),
                settings.CheckboxParameter(name="ShowSingleTent",
                                           label="Show single tent"),
                animate_tent],
            [settings.ValueParameter(
                name="SelectedLevel", label="Level",
                default_value=self._initial_values["SelectedLevel"],
                min_value=0, max_value=self.nlevels-1),
                settings.CheckboxParameter(name="ShowSingleLevel",
                                           label="Show single level"),
                animate_lvl]]

        self.addParameters(
            "Draw Tents",
            settings.SingleChoiceParameter(
                name="Selection",  # label="Draw by",
                options=["all", "by tent number", "by level number"],
                default_value="all",
                sub_parameters=sub_parameters),
            settings.ValueParameter(
                name="ScaleTents", label="Scale",
                default_value=self._initial_values["ScaleTents"],
                min_value=0.0, step=0.2),
            settings.ValueParameter(
                name="ShrinkTents", label="Shrink",
                default_value=self._initial_values["ShrinkTents"],
                min_value=0.0, max_value=1.0, step=0.01))

    @inmain_decorator(True)
    def render(self, settings):
        if not self.active:
            return
        GL.glBindVertexArray(self.vao)

        elements = self.mesh_data.elements[ngsolve.VOL][0]
        prog = getProgram("tents.vert", "tents.geom", "tents.frag",
                          elements=elements, params=settings, scene=self)
        uniforms = prog.uniforms

        uniforms.set('light_ambient', 0.3)
        uniforms.set('light_diffuse', 0.7)

        # set 2d-tent data
        GL.glActiveTexture(GL.GL_TEXTURE2)
        self.tex_tents.bind()
        uniforms.set('tents', 2)

        GL.glActiveTexture(GL.GL_TEXTURE3)
        self.tex_times.bind()
        uniforms.set('times', 3)

        GL.glActiveTexture(GL.GL_TEXTURE5)
        self.tex_tents_color.bind()
        uniforms.set('colors', 5)

        uniforms.set('scale_tents', self.getScaleTents())
        uniforms.set('shrink_tents', self.getShrinkTents())
        if (DrawTents(self.getSelection()) == DrawTents.TENTNR):
            uniforms.set('selected_level', -1)
            uniforms.set('selected_tent', self.getSelectedTent())
            uniforms.set('single_tent', self.getShowSingleTent())
        elif (DrawTents(self.getSelection()) == DrawTents.LEVELNR):
            uniforms.set('selected_level', self.getSelectedLevel())
            uniforms.set('single_level', self.getShowSingleLevel())
            uniforms.set('selected_tent', -1)
        else:
            uniforms.set('selected_level', -1)
            uniforms.set('selected_tent', -1)

        # filled tents
        uniforms.set('wireframe', False)
        GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL)
        GL.glPolygonOffset(0, 0)
        GL.glEnable(GL.GL_POLYGON_OFFSET_FILL)
        GL.glDrawArrays(GL.GL_POINTS, 0, self.lententdata)
        GL.glDisable(GL.GL_POLYGON_OFFSET_FILL)
        GL.glBindVertexArray(0)


class TentSceneSlab(TentScene, settings.ColormapSettings):
    @inmain_decorator(wait_for_return=True)
    def __init__(self, tentslab, **kwargs):
        self._initial_values = {"SelectedTent": 0,
                                "ShowSingleTent": False,
                                "SelectedLevel": 0,
                                "ShowSingleLevel": False,
                                "ShrinkTents": 0.95,
                                "ScaleTents": 1.0}

        mesh = tentslab.mesh
        tentdata, tenttimes, ntents, nlevels = tentslab.DrawPitchedTentsGL()
        self.tentdata = tentdata
        self.tenttimes = tenttimes
        self.nlevels = nlevels
        self.ntents = ntents
        self.lententdata = len(tentdata)//4
        # draw basemesh
        BaseMeshScene.__init__(self, mesh, **kwargs)

    @inmain_decorator(True)
    def update(self):
        BaseMeshScene.update(self)
        # store tent data on GPU memory
        tents = np.array(self.tentdata, dtype=np.int32)
        n = len(tents)
        self.tex_tents.store(tents, GL.GL_INT)
        tenttimes = np.zeros(n, dtype=np.float32) if self.tenttimes is None \
                    else np.array(self.tenttimes, dtype=np.float32)
        self.tex_times.store(tenttimes)
        self.tents_colors = [0, 255, 0, 255]*self.ntents
        self.tex_tents_color.store(self.tents_colors, GL.GL_UNSIGNED_BYTE,
                                   self.ntents)
