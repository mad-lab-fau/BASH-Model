#version 440

in vec3 worldPosition;
in vec3 worldNormal;
in vec3 boneColor;
in float muscleActivity;

uniform vec3 cameraPosition;
uniform bool showNormals;
uniform bool showBodyparts;
uniform bool showMuscleActivity;

out vec4 out_col;

#define BIND_ID_BUFFER_POINTLIGHTS 1
#define BIND_ID_BUFFER_MUSCLEDATA 2

// Point lights
struct PointLight {    
    vec3 position;
    float shininess;
    vec4 ambient;
    vec4 diffuse;
    vec4 specular;
}; 
layout (std430, binding = BIND_ID_BUFFER_POINTLIGHTS) buffer PointLightBuffer {
    PointLight pointLights[];
};



// estimates the normal of a face using the x/y-derivatives
vec3 computeFaceNormal(vec3 wpos) {
	vec3 dFdxPos = dFdx(wpos);
	vec3 dFdyPos = dFdy(wpos);
	vec3 faceNormal = normalize(cross(dFdxPos,dFdyPos));
	return faceNormal;
}


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
	//vec3 normal = world_normal;
	vec3 normal = computeFaceNormal(worldPosition);

	vec3 color = vec3(0.0);
	
	pointLights[0].position = cameraPosition; // fix front light to camera
	for (int i = 0; i < pointLights.length(); i++) {
		color += computeLighting(pointLights[i], worldPosition, normal, cameraPosition);
	}
	
	// muscle activity
	if (showMuscleActivity) {
		if (muscleActivity > 0) {
			color += vec3(muscleActivity, 0.0, 0.0);
			//color += vec3(1.0, 0.0, 0.0);
		}
	}


	// Debug: Normals
	if (showNormals) {
		color = normal;
	}
	// Debug: Show Body Parts
	if (showBodyparts) {
		color = boneColor / 255.0;
	}

	out_col = vec4(color, 1.0);
}
