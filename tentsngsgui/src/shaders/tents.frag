#version 150

{include utils.inc}

uniform bool wireframe;

in VertexData
{
  vec3 pos;
  vec3 normal;
  vec4 color;
  vec3 edgedist;
} inData;

out vec4 FragColor;

void main()
{
  if(wireframe)
    {
      FragColor = vec4(0,0,0,1);
      // float d = min(min(inData.edgedist.x, inData.edgedist.y), inData.edgedist.z);
      // if(d>1e-5) discard;
    }
  else
    {
      FragColor = inData.color;
      
      if (FragColor.a == 0.0)
        discard;
      FragColor.rgb = CalcLight(FragColor.rgb, MV, inData.pos, inData.normal);
    }
}
