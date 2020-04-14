#pragma once
#pragma warning(disable : 4201)

#include "paths.hpp"

#include <exception>
#include <iostream>
#include <string>

#include <atlas/glx/Buffer.hpp>
#include <atlas/glx/Context.hpp>
#include <atlas/glx/ErrorCallback.hpp>
#include <atlas/glx/GLSL.hpp>

#include <fmt/printf.h>
#include <magic_enum.hpp>
#include <atlas/math/Math.hpp>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace atlas;
using Colour = atlas::math::Vector;
using Point = std::array<float, 3>;
static const std::vector<std::string> IncludeDir{ ShaderPath };

class Square;
class Camera;
class PointLight;
class DirectionalLight;

struct OpenGLError : std::runtime_error
{
    OpenGLError(const std::string& what_arg) : std::runtime_error(what_arg) {};
    OpenGLError(const char* what_arg) : std::runtime_error(what_arg) {};
};

class Camera
{
public:
    Camera();
    Camera(glm::vec3 position, glm::vec3 target);

    void calculateRest();
    glm::vec3 mUp;
    glm::vec3 mPosition;
    glm::vec3 mTarget;

protected:
    glm::vec3 mDirection;
    glm::vec3 mFront;
    glm::vec3 mRight;
    glm::vec3 mWorldUp;
};

class PointLight
{
public:
    PointLight();
    PointLight(glm::vec3 point, glm::vec3 diffuse, glm::vec3 pSpecular);

    glm::vec3 mPoint;
    glm::vec3 mDiffuse;
    glm::vec3 mSpecular;
};

class DirectionalLight
{
public:
    DirectionalLight();
    DirectionalLight(glm::vec3 direction, glm::vec3 diffuse, glm::vec3 dSpecular);

    glm::vec3 mDirection;
    glm::vec3 mDiffuse;
    glm::vec3 mSpecular;
private:

};

struct World
{
    Colour background;
    std::vector<Square> scene;
    std::vector<Colour> image;
    PointLight pLight;
    DirectionalLight dLight;
    std::vector<glm::vec3> cubePositions;
    glm::mat4 view;
    Camera cam;
};


class Square
{
public:
    // gl_Position = model * projection * view * transform * vec4(position, 1.0);
    Square();

    void loadShaders();

    void loadDataToGPU(std::array<float, 324> const& vertices);

    void reloadShaders();

    void render(glm::mat4 cameraView, glm::vec3 cubePos);

    void freeGPUData();

    float SpecularStrength;

private:
    // Vertex buffers.
    GLuint mVao;
    GLuint mVbo;
    // Shader data.
    GLuint mVertHandle;
    GLuint mFragHandle;
    GLuint mProgramHandle;
    glx::ShaderFile vertexSource;
    glx::ShaderFile fragmentSource;

};

class Program
{
public:
    Program(int width, int height, std::string title);
    void run();

    void freeGPUData();

private:
    static void errorCallback(int code, char const* message)
    {
        fmt::print("error ({}): {}\n", code, message);
    }

    void createGLContext();

    GLFWwindow* mWindow;
    glx::WindowSettings settings;
    glx::WindowCallbacks callbacks;
};
