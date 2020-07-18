from ipywidgets import DOMWidget, register
import numpy as np
import ngsolve as ngs
import os

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

    def GenerateHTML(self, filename=None):
        import json
        print(" Generating data in Python")
        d = self.GetData()
        print(" Converting data to JSON")
        data = json.dumps(d)
        jscode = "var render_data = {}\n".format(data) + render_js_code
        html = html_template.replace('{render}', jscode)
        if filename is not None:
            open(filename, 'w').write(html)
        return html

    def Draw(self):
        from IPython.display import display
        self.widget = NGSTentsWebGuiWidget()
        # print("self.widget", self.widget)
        d = self.GetData()
        # print("d.keys", d.keys())
        # print("len(d['tent_nrs'])", len(d['tent_nrs']))
        self.widget.value = d
        # print("set self.widget.value")
        # print("self.widget.__dict__", self.widget.__dict__)
        display(self.widget)
        print("called display self.widget")
        return self.widget

    def Redraw(self):
        d = self.GetData(set_minmax=False)
        self.widget.value = d

    def __repr__(self):
        return ""

    def GetData(self):
        d = {}
        lists = self.GetElements2D(d)
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

    def GetElements2D(self, d):
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
        for i, d in enumerate(data):
            nr, layer, vnr, elnr = d
            nrs.append(int(nr))  # convert from int64 to int
            layers.append(int(layer))
            tbots = times[i][:3]
            ttop = times[i][3]
            vertices = []
            eid = ngs.ElementId(ngs.VOL, elnr)
            for j, v in enumerate(self.mesh[eid].vertices):
                vpt = list(self.mesh[v].point)
                # set z-coordinate to (bottom) time of vertex
                vpt.append(tbots[j])
                vertices.append(vpt)
                # central vertex also has a top time
                if v.nr == vnr:
                    vpt_top = vpt[:]
                    vpt_top[2] = ttop
                    vertices.append(vpt_top)
                    center = j
            tentcenters.append(center)
            fcs, nrmls = self.GetFacesAndNormals(vertices)
            faces += fcs
            gvertices += vertices
            normals += nrmls

        return gvertices, normals, faces, tentcenters, nrs, layers

    def GetFacesAndNormals(self, element):
        """
        Given a spacetime element as a list of vertices, compute its
        faces as lists of vertex indices, ordered counter-clockwise from
        outside the element.  Also return the unit normals to each face
        in order corresponding to faces.
        """
        faces = []
        for i in range(4):
            j, k, m = (i+1) % 4, (i+2) % 4, (i+3) % 4
            signs, normals = self.GetSignsAndNormals(element)
            faces.append([j, k, m] if signs[i] > 0 else [j, m, k])
        return faces, normals

    def GetSignsAndNormals(self, element):
        """
        Get the sign determining the orientation of the face vertices
        for each face of the element.
        """
        signs = []
        normals = []
        pts = np.array(element).T  # four column vectors
        vecs = pts - np.roll(pts[:], 1, 1)
        for i in range(4):
            j, k = (i+2) % 4, (i+3) % 4
            normal = np.cross(vecs[:, j], vecs[:, k])
            sign = - np.sign(vecs[:, i].dot(normal))
            signs.append(sign)
            normal = normal / np.linalg.norm(normal) * sign
            normals.append(normal)
        return signs, normals


def Draw(tps, filename='output.html'):
    scene = WebGLScene(tps)
    if _IN_IPYTHON:
        if _IN_GOOGLE_COLAB:
            from IPython.display import display, HTML
            html = scene.GenerateHTML()
            display(HTML(html))
        else:
            print("should render using widgets.DOMWidget")
            # render scene using widgets.DOMWidget
            scene.Draw()
            return scene
    else:
        scene.GenerateHTML(filename=filename)
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
