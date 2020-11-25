#version 440

layout (location = 0) in vec3 in_position;
layout (location = 1) in vec3 in_normal;
layout (location = 2) in vec3 in_color;
layout (location = 3) in ivec4 in_boneID;
layout (location = 4) in vec4 in_weight;

const int NUM_BONES = 20;
uniform mat4 B[NUM_BONES];
uniform mat4 M;
uniform mat4 V;
uniform mat4 P;

out vec3 worldPosition;
out vec3 worldNormal;
out vec3 vertexColor;

void main() {
	mat4 M_bone = B[in_boneID[0]] * in_weight[0] + B[in_boneID[1]] * in_weight[1] + B[in_boneID[2]] * in_weight[2] + B[in_boneID[3]] * in_weight[3];
	
	vec4 T_position = M * M_bone * vec4(in_position, 1.0);
	vec4 T_normal = M * M_bone * vec4(in_normal, 0.0);

	gl_Position = P * V * T_position;
	
	worldPosition = T_position.xyz;
	worldNormal = normalize(T_normal.xyz);
	vertexColor = in_color;
}