#version 330 core

layout (location = 0) in vec3 in_position;

uniform mat4 M;
uniform mat4 V;
uniform mat4 P;

void main() {
	vec4 T_pos = M * vec4(in_position, 1.0);

	gl_Position = P * V * T_pos;
}