import {
  MODULE_NAME, MODULE_VERSION
} from './version';

import * as THREE from 'three';
import dat from 'dat.gui';


function readB64(base64) {
  var binary_string = window.atob(base64);
  var len = binary_string.length;
  var bytes = new Uint8Array(len);
  for (var i = 0; i < len; i++) {
    bytes[i] = binary_string.charCodeAt(i);
  }
  return new Float32Array( bytes.buffer);
}

function setKeys (dst, src) {
  for(var key in dst) {
    if(typeof(dst[key])=="object" && src[key] !== undefined)
      setKeys(dst[key], src[key]);
    else
       dst[key] = src[key];
  }
}

let CameraControls = function(cameraObject, scene, domElement) {
  if ( domElement === undefined ) console.log( 'domElement is undefined' );
  if ( domElement === document )
      console.error('"document" should not be used as the target "domElement".'
       + '  Please use "renderer.domElement" instead.' );
  if ( !cameraObject.isPerspectiveCamera )
      console.error('camera must be perspective camera');

  this.scene = scene;
  this.slab_radius = scene.slab_radius;
  this.center = scene.slab_center.clone();

  this.cameraObject = cameraObject;
  this.pivotObject = scene.pivot;
  this.domElement = domElement;

  this.transmat = new THREE.Matrix4();
  this.rotmat = new THREE.Matrix4();
  this.centermat = new THREE.Matrix4();
  this.transformationmat = new THREE.Matrix4();
  this.scale = 1.0/this.slab_radius;

  this.centermat.makeTranslation(-this.center.x, -this.center.y, -this.center.z);

  this.mode = null;

  this.keys = { LEFT: 37, UP: 38, RIGHT: 39, DOWN: 40, CLOCKWISE: 65,
                COUNTERCLOCKWISE: 83};

  this.rotation_step_degree = 0.05;
  this.pan_step = 0.05;
  this.camera_step = 0.2;

  // not to change from outside
  var changeEvent = { type: 'change' };

  var scope = this;

  this.reset = () => {
    scope.transmat.identity();
    scope.rotmat.identity();
    scope.centermat.identity();
    scope.transformationmat.identity();
    scope.scale = 1.0/this.slab_radius;
    scope.center.copy(scene.slab_center);
    scope.centermat.makeTranslation(
        -this.center.x, -this.center.y, -this.center.z);
    scope.update();
    //scope.scene.setCenterTag();
  }

  this.update = function () {
    var scale_vec = new THREE.Vector3();
    return function update() {
      scale_vec.setScalar(scope.scale);
      scope.pivotObject.matrix.copy(scope.transmat)
           .multiply(scope.rotmat).scale(scale_vec).multiply(scope.centermat);
      const aspect = this.domElement.offsetWidth/this.domElement.offsetHeight;
      this.scene.axes_object.matrixWorld.makeTranslation(-0.85*aspect, -0.85, 0)
          .multiply(scope.rotmat);
      scope.dispatchEvent( changeEvent );
    };
  }()

  this.rotateObject = function () {
    var mat = new THREE.Matrix4();
    return function(axis, rad) {
      mat.makeRotationAxis(axis, rad);
      scope.rotmat.premultiply(mat);
    };
  }();

  this.panObject = function () {
    var mat = new THREE.Matrix4();
    return function(dir, dist) {
      mat.makeTranslation(dist*dir.x, dist*dir.y, dist*dir.z);
      scope.transmat.premultiply(mat);
    };
  }();

  this.updateCenter = function () {
    return function() {
      console.log("set mesh center to", scope.center);
      scope.centermat.makeTranslation(
          -scope.center.x, -scope.center.y, -scope.center.z);
      scope.transmat.identity();
      //scope.scene.setCenterTag();
      scope.update();
    };
  }();

  function keydown(event) {
    var needs_update = false;

    if (event.shiftKey){ // pan
      if (event.keyCode == scope.keys.DOWN) {
        needs_update = true;
        scope.panObject(new THREE.Vector3(0, -1, 0), scope.pan_step)
      } else if (event.keyCode == scope.keys.UP) {
        needs_update = true;
        scope.panObject(new THREE.Vector3(0, 1, 0), scope.pan_step)
      } else if (event.keyCode == scope.keys.LEFT) {
        needs_update = true;
        scope.panObject(new THREE.Vector3(-1, 0, 0), scope.pan_step)
      } else if (event.keyCode == scope.keys.RIGHT) {
        needs_update = true;
        scope.panObject(new THREE.Vector3(1, 0, 0), scope.pan_step)
      }

    } else { // rotate
      if (event.keyCode == scope.keys.DOWN) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(1, 0, 0), scope.rotation_step_degree)
      } else if (event.keyCode == scope.keys.UP) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(-1, 0, 0), scope.rotation_step_degree)
      } else if (event.keyCode == scope.keys.LEFT) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(0, -1, 0), scope.rotation_step_degree)
      } else if (event.keyCode == scope.keys.RIGHT) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(0, 1, 0), scope.rotation_step_degree)
      } else if (event.keyCode == scope.keys.CLOCKWISE) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(0, 0, 1), scope.rotation_step_degree)
      } else if (event.keyCode == scope.keys.COUNTERCLOCKWISE) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(0, 0, -1), scope.rotation_step_degree)
      }
    }

    if(needs_update) {
      event.preventDefault();
      scope.update();
    }

  }

  function onMouseDown(event) {
    if(event.button==0) {
      event.preventDefault();
      scope.mode = "rotate";
    }
    if(event.button==2) {
      event.preventDefault();
      scope.mode = "move";
    }
    event.stopPropagation();
  }

  function onMouseUp(event) {
    scope.mode = null;
    scope.dispatchEvent( changeEvent );
  }


  function onMouseMove(event) {
    var needs_update = false;

    if(scope.mode=="rotate")
      {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(1, 0, 0), 0.01*event.movementY);
        scope.rotateObject(new THREE.Vector3(0, 1, 0), 0.01*event.movementX);
      }

      if(scope.mode=="move")
        {
          needs_update = true;
          scope.panObject(new THREE.Vector3(1, 0, 0), 0.004*event.movementX);
          scope.panObject(new THREE.Vector3(0, -1, 0), 0.004*event.movementY);
        }

        if(needs_update) {
          event.preventDefault();
          scope.update();
        }
  }

  var oldtouch = new THREE.Vector2(0,0);
  var olddelta = 0;
  var touchfirst = true;

  function onTouchStart(event) {
    touchfirst = true;
  }

  function onTouchMove(event) {

    event.preventDefault();

    switch ( event.touches.length ) {
      case 1:
        var pos = new THREE.Vector2(event.touches[0].pageX,
                                    event.touches[0].pageY);
      if (!touchfirst) {
        scope.rotateObject(new THREE.Vector3(1, 0, 0), 0.01*(pos.y-oldtouch.y));
        scope.rotateObject(new THREE.Vector3(0, 1, 0), 0.01*(pos.x-oldtouch.x));
      }
      oldtouch = pos;
      touchfirst = false;
      scope.update();
      break;

      default: // 2 or more
        var dx = event.touches[ 0 ].pageX - event.touches[ 1 ].pageX;
      var dy = event.touches[ 0 ].pageY - event.touches[ 1 ].pageY;
      var delta = Math.sqrt( dx * dx + dy * dy );
      if (!touchfirst) {
        var s = Math.exp(0.01*(delta-olddelta));
        scope.scale *=  s;
      }
      touchfirst = false;
      scope.update();
      olddelta = delta;
      break;
    }
  }

  function wheel(event) {
    event.preventDefault();
    event.stopPropagation();

    let dy = event.deltaY;
    if(event.deltaMode==1) // 1==DOM_DELTA_LINE -> scroll in lines, not pixels
      dy *= 30;

    var s = Math.exp(-0.001*dy);
    scope.scale *=  s ;
    scope.update();
  }

  function contextmenu( event ) {
    event.preventDefault();
  }

  function getPixel(scene, mouse){

  }

  function onDblClick( event ){
    event.preventDefault();
    var rect = scope.domElement.getBoundingClientRect();
    scope.scene.mouse.set(event.clientX-rect.left, event.clientY-rect.top);
    scope.dispatchEvent( changeEvent );

  }

  scope.domElement.addEventListener('dblclick', onDblClick, false);

  window.addEventListener( 'mouseup', onMouseUp, false );
  scope.domElement.addEventListener( 'mousedown', onMouseDown, false );
  scope.domElement.addEventListener( 'contextmenu', contextmenu, false );
  window.addEventListener( 'mousemove', onMouseMove, false );

  scope.domElement.addEventListener( 'touchstart', onTouchStart, false );
  scope.domElement.addEventListener( 'touchmove', onTouchMove, false );

  scope.domElement.addEventListener( 'wheel', wheel, false );


  if ( scope.domElement.tabIndex === - 1 ) {

    scope.domElement.tabIndex = 0;

  }

  this.reset();
}; // end CameraControls function definition



export class Scene {
  render_data: any;
  scene: THREE.Scene;
  renderer: THREE.WebGLRenderer;
  camera: THREE.PerspectiveCamera;
  ortho_camera: THREE.OrthographicCamera;
  light: any;
  amblight: any;
  container: any;
  stats: any;

  gui: any;
  gui_status_default: any;
  gui_status: any;
  gui_functions: any;
  axes_object: any;

  context: any;
  trafo: any;
  render_target: any;
  mouse: THREE.Vector2;

  last_frame_time: number;

  tent_colors: any;
  tent_elements: any;
  slab_center: THREE.Vector3;
  slab_radius: number;
  ntents: number;
  nlayers: number;

  requestId: number;
  have_webgl2: boolean;

  pivot: THREE.Group;

  label_style: string;

  controls: any;
  element: any;

  version_object: any;

  constructor() {
    this.gui_status_default = {
      Light: { ambient: 2, directional: .5},
      Material: { shininess: 30},
      Misc: { stats: "-1", reduce_subdivision: false,
              "version": true, "axes": true, "colormap": true },
      display_by: "layers",
      z_scaling: 4.0,
      shrink_tents: 0.98,
      tent: 0,
      layer: 0,
    };
    // deep-copy settings
    this.gui_status = JSON.parse(JSON.stringify(this.gui_status_default));
    this.gui_functions = { };

    this.tent_elements = [];
    this.tent_colors = [0xe6194b, 0x3cb44b, 0xffe119, 0x4363d8, 0xf58231,
                   0x911eb4, 0x46f0f0, 0xf032e6, 0xbcf60c, 0xfabebe,
                   0x008080, 0xe6beff, 0x9a6324, 0xfffac8, 0x800000,
                   0xaaffc3, 0x808000, 0x000075, 0x808080, 0x000000];


    this.have_webgl2 = false;

    this.label_style  = '-moz-user-select: none; -webkit-user-select: none;'
                        +' -ms-user-select:none; onselectstart="return false;';
    this.label_style += 'onmousedown="return false; user-select:none;'
                        +' -o-user-select:none;unselectable="on";';
    this.label_style += 'position: absolute; z-index: 100; display:block;';
    this.requestId = 0;
  }

  setGuiSettings (settings) {
    console.log("in gui settings");
    setKeys(this.gui_status, settings);
    for (var i in this.gui.__controllers)
      this.gui.__controllers[i].updateDisplay();
    for (var f in this.gui.__folders) {
      const folder = this.gui.__folders[f];
      for (var i in folder.__controllers)
        folder.__controllers[i].updateDisplay();
    }
    this.animate();
  }

  onResize() {
    const w = this.element.parentNode.clientWidth;
    const h = this.element.parentNode.clientHeight;

    const aspect = w/h;
    this.ortho_camera = new THREE.OrthographicCamera( -aspect, aspect,
                                                      1.0, -1.0, -100, 100 );
    this.camera.aspect = aspect;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize( w, h );
    if (this.controls)
    {
      this.controls.update();
    }
    this.animate();
  }

  init (element, render_data)
  {
    this.last_frame_time = new Date().getTime();
    this.render_data = render_data;
    this.element = element;
    console.log("THREE", THREE);
    console.log("dat", dat);
    // console.log("Stats", Stats);

    CameraControls.prototype = Object.create( THREE.EventDispatcher.prototype );
    CameraControls.prototype.constructor = CameraControls;

    this.ntents = render_data.ntents;
    this.nlayers = render_data.nlayers;
    this.slab_radius = render_data.slab_radius;
    this.slab_center = new THREE.Vector3().fromArray(render_data.slab_center);

    var canvas = document.createElement( 'canvas' );

    var gl2 = canvas.getContext('webgl2');

    if (gl2) {
      console.log('webgl2 is supported!');
      this.context = canvas.getContext( 'webgl2', { alpha: false } );
      this.have_webgl2 = true;
    }
    else
    {
      console.log('your browser/OS/drivers do not support WebGL2');
      this.context = canvas.getContext( 'webgl', { alpha: false } );
    }

    this.renderer = new THREE.WebGLRenderer( { canvas: canvas,
                                               context: this.context } );
    this.renderer.autoClear = false;
    console.log("Renderer", this.renderer);

    this.render_target = new THREE.WebGLRenderTarget( window.innerWidth,
                                                      window.innerHeight );
    this.render_target.texture.format = THREE.RGBAFormat;
    this.render_target.texture.type = THREE.FloatType;

    //this is to get the correct pixel detail on portable devices
    this.renderer.setPixelRatio( window.devicePixelRatio );

    //and this sets the canvas' size.
    this.renderer.setSize( this.element.offsetWidth, this.element.offsetHeight );
    this.renderer.setClearColor( 0xffffff, 1 );

    this.container = document.createElement( 'div' );
    element.appendChild( this.container );

    this.container.appendChild( this.renderer.domElement );

    // label with NGSolve version at right lower corner
    this.version_object = document.createElement("div");
    var style = 'bottom: 10px; right: 10px';
    this.version_object.setAttribute("style",this.label_style+style);
    var version_text = document.createTextNode(
        "NGSolve " + render_data.ngsolve_version);
    this.version_object.appendChild(version_text)
    this.container.appendChild(this.version_object);


    this.scene = new THREE.Scene();
    // added to fix background issue on old Mac OS
    this.scene.background = new THREE.Color(0xffffff);
    console.log("scene", this.scene);
    this.axes_object = new THREE.AxesHelper(0.15);
    this.axes_object.matrixAutoUpdate = false;

    this.pivot = new THREE.Group();
    this.pivot.matrixAutoUpdate = false;

    this.camera = new THREE.PerspectiveCamera(
      40,                                         //FOV
      this.element.offsetWidth/ this.element.offsetHeight, // aspect
      1,                                          //near clipping plane
      100                                         //far clipping plane
    );

    this.camera.position.set( 0.0, 0.0, 3 );

    window.addEventListener( 'resize', ()=>this.onResize(), false );

    var color = 0xFFFFFF;
    var intensity = 0.5;
    this.light = new THREE.DirectionalLight(color, intensity);
    this.light.position.set(0,1,3);
    this.scene.add(this.light);

    var ambintensity = 2;
    this.amblight = new THREE.AmbientLight( 0x404040, ambintensity ); // soft white light
    this.scene.add( this.amblight );

    this.trafo = new THREE.Vector2(
        1.0/2.0/(this.slab_center.length()+this.slab_radius), 1.0/2.0);

    this.mouse = new THREE.Vector2(0.0, 0.0);

    let gui = new dat.GUI({autoplace: false, closeOnTop: true});
    let gui_container = document.createElement( 'div' );
    gui_container.setAttribute(
        "style",
        'position: absolute; z-index: 200; display:block; right: 0px; top: 0px');
    gui_container.appendChild(gui.domElement);
    this.container.appendChild(gui_container);

    this.gui = gui;
    console.log("GUI", gui);
    let gui_status = this.gui_status;
    console.log("gui_status", gui_status);
    let animate = ()=>this.animate();

    gui.add(gui_status, "display_by", ["tents", "layers"])
       .name("display by").onChange(animate);
    gui.add(gui_status, "tent", 0, this.ntents-1).step(1).onChange(animate);
    gui.add(gui_status, "layer", 0, this.nlayers-1).step(1).onChange(animate);
    gui.add(gui_status, "z_scaling", 0.01, 100)
       .name("z-scaling").onChange(animate);
    gui.add(gui_status, "shrink_tents", 0.8, 1.0).step(0.0001)
       .name("shrink tents").onChange(animate);

    let gui_light = gui.addFolder("Light");
    gui_light.add(gui_status.Light, "ambient", 0.0, 3.0).onChange(animate);
    gui_light.add(gui_status.Light, "directional", 0.0, 1.0).onChange(animate);

    let gui_material = gui.addFolder("Material");
    gui_material.add(gui_status.Material, "shininess", 0, 100).onChange(animate);

    let gui_misc = gui.addFolder("Misc");
    let gui_functions = this.gui_functions;
    gui_functions['reset settings'] = () =>{
      this.setGuiSettings(this.gui_status_default);
    };
    gui_functions['store settings'] = () => {
      document.cookie = "gui_status="+btoa(JSON.stringify(gui_status)) +
          ";SameSite=Lax";
    };
    gui_functions['load settings'] = () =>{
      var name = "gui_status="
      var decodedCookie = decodeURIComponent(document.cookie);
      var ca = decodedCookie.split(';');
      for(var i = 0; i <ca.length; i++) {
        var c = ca[i];
        c.trim();
        if (c.indexOf(name) == 0) {
          const s = JSON.parse(atob(c.substring(name.length, c.length)));
          this.setGuiSettings(s);
        }
      }
    };
    gui_misc.add(gui_functions, "reset settings");
    gui_misc.add(gui_functions, "store settings");
    gui_misc.add(gui_functions, "load settings");

    gui_misc.add(gui_status.Misc, "axes").onChange(animate);
    gui_misc.add(gui_status.Misc, "version").onChange(value => {
      this.version_object.style.visibility = value ? "visible" : "hidden";
    });

    gui_functions['fullscreen'] = () =>{
      let elem = this.element.parentNode;

      if (elem.requestFullscreen) {
        elem.requestFullscreen();
      } else if(elem.webkitRequestFullScreen) {
        // Webkit (works in Safari and Chrome Canary)
        elem.webkitRequestFullScreen();
      }else if(elem.mozRequestFullScreen) {
        // Firefox
        elem.mozRequestFullScreen();
      }
    };
    gui.add(gui_functions, "fullscreen");

    gui_functions['reset'] = ()=> {
      this.controls.reset();
    };
    gui.add(gui_functions, "reset").onChange(animate);

    gui_functions['update center'] = ()=> {
      this.controls.updateCenter();
    };
    gui.add(gui_functions, "update center").onChange(animate);

    this.scene.add( this.pivot );

    this.controls = new CameraControls(this.camera, this,
                                       this.renderer.domElement );

    console.log(this.controls);
    this.controls.addEventListener('change', animate);

    this.updateRenderData(render_data);
    setTimeout(()=> this.onResize(), 0);
  }

  // called on scene.Redraw() from Python
  updateRenderData(render_data)
  {
    this.render_data = render_data;
    this.setRenderData(render_data);
  }

  setRenderData(render_data)
  {
    let faces = render_data.faces
    let nrs = render_data.tent_nrs
    let layers = render_data.tent_layers
    let vs = readB64(render_data.tent_el_vertices);
    let tent_centers = render_data.tent_centers;
    let nmls = readB64(render_data.face_normals);

    let vi = 0; // global vertex index
    let fi = 0; // global face index
    let ni = 0; // global normals index
    for (let e=0; e < nrs.length; e++)
    {
      let verts = [];       // vertices for each face of element
      let vertnmls = [];    // face normal repeated for each vertex
      // we use the tent center with zero z-coordinate as the base point
      // for each tent element to allow shrinking and scaling without
      // unwanted translation.
      let ctr = tent_centers[e];
      let vctr0 = vs[(vi+ctr)*3];
      let vctr1 = vs[(vi+ctr)*3+1];
      // get the 4 vertices for the element
      let pos = []
      for (let v=0; v < 4; v++)
      {
         pos[v] = vs.slice((vi+v)*3, (vi+v+1)*3);
         pos[v][0] = pos[v][0]-vctr0;
         pos[v][1] = pos[v][1]-vctr1;
      }
      vi += 4
      // for each face, get its normal, then append the three vertices
      // indexed by the face and the normal (3 times)
      for (let f=0; f < 4; f++)
      {
        let face = faces[fi+f];
        let nml = nmls.slice((ni+f)*3, (ni+f+1)*3);
        // push the face vertices in the correct order, then the normal
        for (let v=0; v< 3; v++)
        {
          let vert = pos[face[v]];
          verts.push(...vert);
          vertnmls.push(...nml);
        }
      }
      fi += 4;
      ni += 4;
      let tentGeom = new THREE.BufferGeometry();
      tentGeom.setAttribute(
         'position', new THREE.BufferAttribute(new Float32Array(verts), 3));
      tentGeom.setAttribute(
         'normal', new THREE.BufferAttribute(new Float32Array(vertnmls), 3));

      let nr = nrs[e];
      let layer = layers[e];
      const material = new THREE.MeshPhongMaterial({
        color: this.tent_colors[nr % this.tent_colors.length],
      });
      var elem = new THREE.Mesh(tentGeom, material);
      elem.translateX(vctr0).translateY(vctr1)
      elem.userData = {nr :nr, layer: layer};
      this.tent_elements.push(elem);
      this.pivot.add(elem);
    }
    this.animate();
  }


  animate () {
    // Don't request a frame if another one is currently in the pipeline
    if(this.requestId === 0)
      this.requestId = requestAnimationFrame( ()=>this.render() );

    //   stats.update();
  }

  render() {
    let now = new Date().getTime();
    let frame_time = 0.001*(new Date().getTime() - this.last_frame_time );

    this.requestId = 0;

    if(this.ortho_camera === undefined)
      return; // not fully initialized yet

    let gst = this.gui_status;
    let colors = this.tent_colors;
    this.tent_elements.forEach(function(el) {
      let nr = el.userData.nr;
      let layer = el.userData.layer;
      let num = (gst.display_by === "tents" ? nr : layer);
      el.material.color.setHex( colors[num % colors.length] );
      el.material.shininess = gst.Material.shininess;

      el.visible = (gst.display_by === "tents" && nr <= gst.tent) ||
                     (gst.display_by === "layers" && layer <= gst.layer);
      el.scale.x = gst.shrink_tents;
      el.scale.y = gst.shrink_tents;
      el.scale.z = gst.z_scaling * gst.shrink_tents;
    });

    this.light.intensity = gst.Light.directional;
    this.amblight.intensity = gst.Light.ambient;

    this.axes_object.visible = gst.Misc.axes;

    this.renderer.render( this.scene, this.camera );

    if(this.axes_object && gst.Misc.axes)
      this.renderer.render( this.axes_object, this.ortho_camera );

    this.last_frame_time = now;
  }

}
