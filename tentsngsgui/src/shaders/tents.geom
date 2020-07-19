#version 150 

{include utils.inc}

uniform sampler1D colors;
uniform int selected_level;
uniform bool single_level;
uniform int selected_tent;
uniform bool single_tent;
uniform float scale_tents;
uniform float shrink_tents;
uniform isamplerBuffer tents;
uniform samplerBuffer times;
    
layout(points) in;
layout(triangle_strip, max_vertices=12) out;

in VertexData
{
  flat int element;
} inData[];

out VertexData
{
  vec3 pos;
  vec3 normal;
  vec4 color;
  vec3 edgedist;
} outData;

struct Tent
{
  vec3 pos[4];
  vec3 center;

  int vertex;
  int tentnr;
  int level;
  int elnr;

  bool is_visible;
};

Tent getTent2d(Mesh mesh, int ei ) {
  Tent el;
  
  el.tentnr = texelFetch(tents, ei).r;
  el.level = texelFetch(tents, ei).g;
  el.vertex = texelFetch(tents, ei).b;
  el.elnr = texelFetch(tents, ei).a;

  ELEMENT_TYPE trig = getElement(el.elnr);
  
  vec4 tbot = texelFetch(times, ei);
  // top time of central vertex is last element in times
  float ttop = tbot.a;
  el.pos[3] = texelFetch(mesh.vertices, el.vertex).xyz;
  el.pos[3].z = scale_tents*ttop;

  int offset = mesh.offset + ELEMENT_SIZE*el.elnr;
  for (int i=0; i<3; i++) {
    el.pos[i] = trig.pos[i];
    el.pos[i].z = scale_tents*tbot[i];

    int v = texelFetch(mesh.elements, offset+i+2).r;
    if(v == el.vertex)
    	 el.center = el.pos[i];
  }

  for (int i=0; i<4; i++) {
    el.pos[i] = shrink_tents*el.pos[i] + (1.0-shrink_tents)*el.center;
  }
  el.is_visible = CalcClipping(el.center);
  return el;
}

Tent getTent3d(Mesh mesh, int ei ) {
  Tent el;

  el.tentnr = texelFetch(tents, ei).r;
  el.level = texelFetch(tents, ei).g;
  el.vertex = texelFetch(tents, ei).b;
  el.elnr = texelFetch(tents, ei).a;

  ELEMENT_TYPE tet = getElement(el.elnr);
  int offset = mesh.offset + ELEMENT_SIZE*el.elnr;
  for (int i=0; i<4; i++) {
    el.pos[i] = tet.pos[i];
    int v = texelFetch(mesh.elements, offset+i+2).r;
    if(v == el.vertex)
      el.center = el.pos[i];
  }
  for (int i=0; i<4; i++) {
    el.pos[i] = shrink_tents*el.pos[i] + (1.0-shrink_tents)*el.center;
  }
  el.is_visible = CalcClipping(el.center);
  return el;
}

void AddPoint( vec3 pos ) {
  outData.pos = pos;
  gl_Position = P * MV * vec4(outData.pos, 1);
  EmitVertex();
}

void DrawTrig( Tent el, int face_index ) {
  vec3 pos[3];
  int cnt=0; 
  for( int i=0; i<4; i++ )
    if(i!=face_index)
      {
        pos[cnt] = el.pos[i];
        cnt++;
      }
  vec3 normal = cross(pos[1]-pos[0], pos[2]-pos[0]);
  vec3 center = 0.25*(el.pos[0]+el.pos[1]+el.pos[2]+el.pos[3]);
  if(dot(pos[0]-center, normal)<0)
    outData.normal = (-1)*normal;
  else
    outData.normal = normal;

  float val = 1.6180339887;
  if( selected_level<0 )
    val *= el.tentnr;
  else
    val *= el.level;
  outData.color = vec4(MapColor(val - int(val)),1);
  AddPoint(pos[0]);
  AddPoint(pos[1]);
  AddPoint(pos[2]);
  EndPrimitive();
}

void DrawTent(Tent el) {
  bool draw = false;
  if( selected_tent<0 && selected_level<0 )
    draw = true;
  else if(selected_level<0)
    {
      if(selected_tent == el.tentnr || (el.tentnr < selected_tent && !single_tent))
        draw = true;
    }
  else if(selected_tent<0)
    {
      if(selected_level == el.level || (el.level < selected_level && !single_level))
        draw = true;
    }

  if(draw && el.is_visible )
    {
      DrawTrig(el,0);
      DrawTrig(el,1);
      DrawTrig(el,2);
      DrawTrig(el,3);
    }
}

void main() {
  
  if(mesh.dim==2) {
    Tent el = getTent2d(mesh, inData[0].element);
    DrawTent(el);
  }
  if(mesh.dim==3) {
    Tent el = getTent3d(mesh, inData[0].element);
    DrawTent(el);
  }
}
