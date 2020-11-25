#version 440

layout (location = 0) in vec3 in_position;
layout (location = 1) in vec3 in_normal;
layout (location = 2) in vec3 in_boneColor;
layout (location = 3) in float in_muscleActivity;

uniform mat4 M;
uniform mat4 V;
uniform mat4 P;

out vec3 worldPosition;
out vec3 worldNormal;
out vec3 boneColor;
out float muscleActivity;

void main() {
	vec4 T_position = M * vec4(in_position, 1.0);
	vec4 T_normal = M * vec4(in_normal, 0.0);

	gl_Position = P * V * T_position;
	
	worldPosition = T_position.xyz;
	worldNormal = normalize(T_normal.xyz);
	boneColor = in_boneColor;
	muscleActivity = in_muscleActivity;
}