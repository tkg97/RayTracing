#version 330 core

// Ouput data
out vec3 color;

uniform vec3 in_color;

void main()
{ 
	color = in_color;
}