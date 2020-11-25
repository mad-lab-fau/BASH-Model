#version 330 core

uniform vec3 color;

out vec4 out_col;

void main (void) {
	out_col = vec4(color, 1.0);
}