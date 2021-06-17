from time import time
import os

from ipywidgets import DOMWidget, register
import ngsolve as ngs
import numpy as np

# the build script fills the contents of the variables below
render_js_code = ""
widgets_version = ""


def _jupyter_nbextension_paths():
    return [
        {
            "section": "notebook",
            "src": "nbextension/static",
            "dest": "ngstents_jupyter_widgets",
            "require": "ngstents_jupyter_widgets/extension",
        }
    ]


try:
    __IPYTHON__
    _IN_IPYTHON = True
except NameError:
    _IN_IPYTHON = False

try:
    import google.colab
    _IN_GOOGLE_COLAB = True
except ImportError:
    _IN_GOOGLE_COLAB = False

# <script src="https://cdn.jsdelivr.net/npm/
#              three@0.115.0/build/three.min.js"></script>
# <script src="https://cdnjs.cloudflare.com
#              /ajax/libs/dat-gui/0.7.7/dat.gui.js"></script>
# <script src="https://cdnjs.cloudflare.com
#              /ajax/libs/stats.js/r16/Stats.min.js"></script>

html_template = """
<!DOCTYPE html>
<html>
  <head>
      <meta content="text/html;charset=utf-8" http-equiv="Content-Type">
      <meta content="utf-8" http-equiv="encoding">
      <title>NGS-Tents Visualization</title>
      <meta name='viewport' content='width=device-width, user-scalable=no'/>
      <style>
          body{
                margin:0;
                overflow:hidden;
          }
          canvas{
                cursor:grab;
                cursor:-webkit-grab;
                cursor:-moz-grab;
          }
          canvas:active{
                cursor:grabbing;
                cursor:-webkit-grabbing;
                cursor:-moz-grabbing;
          }
      </style>
  </head>
  <body>
    <script src="https://requirejs.org/docs/release/2.3.6/minified/require.js">
    </script>
    <script>
          {render}

          require(["ngstents_jupyter_widgets"], ngs=>
          {
              let scene = new ngs.Scene();
              scene.init(document.body, render_data);
          });
    </script>
  </body>
</html>
"""


class WebGLScene:
    def __init__(self, tps):
        self.mesh = tps.mesh
        self.tps = tps

    def GenerateHTML(self, filename=None, show_timing=False):
        import json
        print(" Generating data in Python")
        d = self.GetData(show_timing=show_timing)
        print(" Converting data to JSON")
        data = json.dumps(d)
        jscode = "var render_data = {}\n".format(data) + render_js_code
        html = html_template.replace('{render}', jscode)
        if filename is not None:
            open(filename, 'w').write(html)
        return html

    def Draw(self, show_timing=False):
        from IPython.display import display
        self.widget = NGSTentsWebGuiWidget()
        d = self.GetData(show_timing=show_timing)
        self.widget.value = d
        display(self.widget)
        return self.widget

    def Redraw(self):
        d = self.GetData(set_minmax=False)
        self.widget.value = d

    def __repr__(self):
        return ""

    def GetData(self, show_timing=False):
        d = {}
        lists = self.GetElements2D(d, show_timing=show_timing)
        vertices, normals, faces, tentcenters, nrs, layers = lists
        vertices = np.array(vertices)
        normals = np.array(normals)
        vmax, vmin = vertices.max(axis=0), vertices.min(axis=0)

        d['ngsolve_version'] = ngs.__version__
        # distinct element vertices
        d['tent_el_vertices'] = encodeData(vertices)
        # normals for each face
        d['face_normals'] = encodeData(normals)
        # local (element level) indices of vertices for each face CCW
        d['faces'] = faces
        # local central vertices of tents
        d['tent_centers'] = tentcenters
        d['tent_nrs'] = nrs
        d['tent_layers'] = layers
        d['slab_center'] = list((vmin+vmax)/2)
        d['slab_radius'] = np.linalg.norm(vmax-vmin)/2
        return d

    def GetElements2D(self, d, show_timing=False):
        """
        Return six lists
        1. vertices: the distinct vertices for each element, where
            each vertex is a list of 3 coordinates.
        2. normals: the 4 face normal vectors for each element, where
            each vector is a list of 3 coordinates.
        3. faces: indices defining the vertices for each
            face of each element in CCW order
        4. tentcenters: each element is the index of the central vertex
            in the vertices list for the element.  This is used to set a
            base point for each tent in threejs so that shrinking and
            scaling do not cause inadvertent translation of the element.
        5. nrs: tent number for each element
        6. layers: layer for each element

        We can cache the spatial data if there are performance issues,
        but doing it this way ensures that the times match the vertices.
        """
        mesh = self.mesh
        msdict = dict((item[0] for item in
                       mesh.GetPeriodicNodePairs(ngs.VERTEX)))
        data, times, ntents, nlayers = self.tps.DrawPitchedTentsGL()
        d['ntents'] = ntents
        d['nlayers'] = nlayers
        data = np.array(data).reshape(-1, 4)
        times = np.array(times).reshape(-1, 4)
        gvertices = []
        tentcenters = []
        nrs = []
        layers = []
        faces = []
        normals = []
        setup_times = []
        tent_times = []
        face_times = []
        for i, d in enumerate(data):
            start = time()
            nr, layer, vnr, elnr = d
            nrs.append(int(nr))  # convert from int64 to int
            layers.append(int(layer))
            tbots = times[i][:3]
            ttop = times[i][3]
            vertices = []
            eid = ngs.ElementId(ngs.VOL, elnr)
            setup_times.append(time()-start)
            start = time()
            for j, v in enumerate(mesh[eid].vertices):
                vpt = list(mesh[v].point)
                # set z-coordinate to (bottom) time of vertex
                vpt.append(tbots[j])
                vertices.append(vpt)
                # central vertex also has a top time
                # handle case when vnr is master vertex and not in element
                if v.nr == vnr or (vnr in msdict and v.nr == msdict[vnr]):
                    vpt_top = vpt[:]
                    vpt_top[2] = ttop
                    vertices.append(vpt_top)
                    center = j
            tentcenters.append(center)
            tent_times.append(time()-start)
            start = time()
            fcs, nrmls = self.GetFacesAndNormals(vertices)

            faces.extend(fcs)
            gvertices.extend(vertices)
            normals.extend(nrmls)
            face_times.append(time()-start)
        if show_timing:
            print("setup: {:.6f}, tents: {:.6f} faces: {:.6f}".format(
                sum(setup_times), sum(tent_times), sum(face_times)))

        return gvertices, normals, faces, tentcenters, nrs, layers

    def GetFacesAndNormals(self, element):
        """
        Compute the outward normal vectors for each face of the element
        and return the normal vectors along with the faces, where a 'face'
        is a list of vertex indices ordered counterclockwise from outside.
        """
        pts = np.array(element)  # four row vectors of length 3
        vecs = pts - np.roll(pts[:], 1, 0)
        normals = np.cross(np.roll(vecs[:], 2, 0), np.roll(vecs[:], 1, 0))
        signs = - np.sign((vecs*normals).sum(1))
        normals /= ((normals**2).sum(1)**(1./2)*signs).reshape(4, 1)
        faces = [[1, 2, 3], [2, 3, 0], [3, 0, 1], [0, 1, 2]]
        for i, s in enumerate(signs):
            if s < 0:
                faces[i].reverse()
        return faces, normals.tolist()


def Draw(tps, filename='output.html', show_timing=False):
    scene = WebGLScene(tps)
    if _IN_IPYTHON:
        if _IN_GOOGLE_COLAB:
            from IPython.display import display, HTML
            html = scene.GenerateHTML(show_timing=show_timing)
            display(HTML(html))
        else:
            # render scene using widgets.DOMWidget
            scene.Draw(show_timing=show_timing)
            return scene
    else:
        scene.GenerateHTML(filename=filename, show_timing=show_timing)
        return scene


@register
class NGSTentsWebGuiWidget(DOMWidget):
    from traitlets import Dict, Unicode
    _view_name = Unicode('NGSTentsView').tag(sync=True)
    _view_module = Unicode('ngstents_jupyter_widgets').tag(sync=True)
    _view_module_version = Unicode(widgets_version).tag(sync=True)
    value = Dict({"ngstents_version": '0.0.0'}).tag(sync=True)


def encodeData(array):
    from base64 import b64encode
    values = np.array(array.flatten(), dtype=np.float32)
    res = b64encode(values).decode("ascii")
    return res


_jupyter_lab_extension_path = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "labextension")


def howtoInstallJupyterLabextension():
    print("""# To install jupyter lab extension:
jupyter labextension install --clean {labdir}
""".format(labdir=_jupyter_lab_extension_path))
