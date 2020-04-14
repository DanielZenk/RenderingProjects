#version 450 core

struct DirLight {
    vec3 direction;
	
    vec3 diffuse;
    vec3 specular;
};


struct PointLight {
	vec3 diffuse;
	vec3 specular;
};

in vec3 vertexColour;
in vec3 Normal;
in vec3 FragPos;
in vec3 LightPos;
out vec4 fragColour;

uniform DirLight dirLight;
uniform PointLight pointLight;
uniform vec3 cameraPos;
uniform float specularStr;

vec3 CalcDirLight(DirLight light, vec3 normal, vec3 viewDir, float specularStr);
vec3 CalcPointLight(PointLight light, vec3 normal, vec3 fragPos, vec3 viewDir, float specularStr, vec3 lightPos);


void main() {
	vec3 norm = normalize(Normal);
	vec3 lightDir = normalize(LightPos - FragPos); 
	
	float diff = max(dot(norm, lightDir), 0);
	vec3 diffuse = diff * vec3(1.0, 1.0, 1.0);

	vec3 viewDir = normalize(cameraPos - FragPos);
	vec3 result = CalcPointLight(pointLight, norm, FragPos, viewDir, specularStr, LightPos);

	result += CalcDirLight(dirLight, norm, viewDir, specularStr);
    result *= vertexColour;
	fragColour = vec4(result, 1.0);
}

vec3 CalcDirLight(DirLight light, vec3 normal, vec3 viewDir, float specularStr)
{
    vec3 lightDir = normalize(-light.direction);
    // diffuse shading
    float diff = max(dot(normal, lightDir), 0.0);
    // specular shading
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 diffuse = light.diffuse * diff * vec3(1.0);
    vec3 specular = specularStr * spec * vec3(1.0);
    return (diffuse + specular);
}

vec3 CalcPointLight(PointLight light, vec3 normal, vec3 fragPos, vec3 viewDir, float specularStr, vec3 lightPos)
{
    vec3 lightDir = normalize(lightPos - fragPos);
    // diffuse shading
    float diff = max(dot(normal, lightDir), 0.0);
    // specular shading
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    
    // combine results
    vec3 diffuse = light.diffuse * diff * vec3(1.0);
    vec3 specular = specularStr * spec * vec3(1.0);
    return (diffuse + specular);
}