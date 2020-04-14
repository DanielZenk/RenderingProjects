#version 450 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 colour;
layout(location = 2) in vec3 normal;

out vec3 vertexColour;
out vec3 Normal;
out vec3 FragPos;
out float Radiance;
out vec3 LightPos;

uniform mat4 transform;
uniform vec3 lightPos;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 model;
uniform float radiance;

void main() {
	FragPos = vec3(model * transform * vec4(position, 1.0));
	gl_Position = projection * view  * model * transform * vec4(position, 1.0);
	vertexColour = colour;
	LightPos = vec3(model * view *  vec4(lightPos, 0.0));
	Normal = vec3(transform * vec4(normal,1));
	Radiance = radiance;
}