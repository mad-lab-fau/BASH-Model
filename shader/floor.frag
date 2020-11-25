#version 440

in vec3 worldPosition;
in vec3 worldNormal;

uniform vec3 cameraPosition;

out vec4 out_col;

#define BIND_ID_LIGHTBUFFER 1

// Point lights
struct PointLight {    
    vec3 position;
    float shininess;
    vec4 ambient;
    vec4 diffuse;
    vec4 specular;
}; 
layout (std430, binding = BIND_ID_LIGHTBUFFER) buffer LightBuffer {
    PointLight pointLights[];
};


// phong-blinn lighting in world space
vec3 computeLighting(PointLight light, vec3 wpos, vec3 normal, vec3 campos) {
	vec3 lightDir = normalize(light.position - wpos); 
	vec3 viewDir = normalize(campos - wpos);
	vec3 reflection = normalize(-reflect(lightDir, normal));
   
	vec3 ambient = light.ambient.rgb;
	vec3 diffuse = clamp(light.diffuse.rgb * max(dot(normal, lightDir), 0.0), 0.0, 1.0);
	vec3 specular = clamp(light.specular.rgb * pow(max(dot(reflection, viewDir), 0.0), light.shininess), 0.0, 1.0);

	return (ambient + diffuse + specular);
} 

void main (void) {
	vec3 normal = worldNormal;

	vec3 color = vec3(0.0);
	for (int i = 0; i < pointLights.length(); i++) {
		color += computeLighting(pointLights[i], worldPosition, normal, cameraPosition);
	}


	out_col = vec4(color, 1.0);
}
