#include "main.hpp"


#define SCR_HEIGHT 720
#define SCR_WIDTH 1280


std::shared_ptr<World> world{ std::make_shared<World>() };
bool directional = true;
bool specular = true;

// ===---------------SQUARE-----------------===
Square::Square()
{
    mProgramHandle = glCreateProgram();
    mVertHandle = glCreateShader(GL_VERTEX_SHADER);
    mFragHandle = glCreateShader(GL_FRAGMENT_SHADER);
}

void Square::loadShaders()
{
    std::string shaderRoot{ ShaderPath };
    vertexSource = glx::readShaderSource(shaderRoot + "Square.vert", IncludeDir);
    fragmentSource = glx::readShaderSource(shaderRoot + "Square.frag", IncludeDir);

    if (auto result{ glx::compileShader(vertexSource.sourceString, mVertHandle) };
        result) {
        throw OpenGLError(*result);
    }

    if (auto result{ glx::compileShader(fragmentSource.sourceString, mFragHandle) };
        result) {
        throw OpenGLError(*result);
    }

    glAttachShader(mProgramHandle, mVertHandle);
    glAttachShader(mProgramHandle, mFragHandle);

    if (auto result = glx::linkShaders(mProgramHandle); result) {
        throw OpenGLError(*result);
    }
}

void Square::loadDataToGPU(
    [[maybe_unused]] std::array<float, 324> const& vertices)
{
    glCreateBuffers(1, &mVbo);

    glNamedBufferStorage(
        mVbo, glx::size<float>(vertices.size()), vertices.data(), 0
    );
    glCreateVertexArrays(1, &mVao);


    glEnableVertexArrayAttrib(mVao, 0);
    glEnableVertexArrayAttrib(mVao, 1);
    glEnableVertexArrayAttrib(mVao, 2);

    glVertexArrayVertexBuffer(mVao, 0, mVbo, 0, glx::stride<float>(9));

    glVertexArrayAttribFormat(
        mVao, 0, 3, GL_FLOAT, GL_FALSE, glx::relativeOffset<float>(0)
    );

    glVertexArrayAttribFormat(
        mVao, 1, 3, GL_FLOAT, GL_FALSE, glx::relativeOffset<float>(3)
    );

    glVertexArrayAttribFormat(
        mVao, 2, 3, GL_FLOAT, GL_FALSE, glx::relativeOffset<float>(6)
    );

    glVertexArrayAttribBinding(mVao, 0, 0);
    glVertexArrayAttribBinding(mVao, 1, 0);
    glVertexArrayAttribBinding(mVao, 2, 0);
}

void Square::reloadShaders()
{
    if (glx::shouldShaderBeReloaded(vertexSource)) {
        glx::reloadShader(mProgramHandle, mVertHandle, vertexSource, IncludeDir);
    }

    if (glx::shouldShaderBeReloaded(fragmentSource)) {
        glx::reloadShader(mProgramHandle, mFragHandle, fragmentSource, IncludeDir);
    }
}


void Square::render(glm::mat4 cameraView, [[maybe_unused]]glm::vec3 cubePos)
{

    glUseProgram(mProgramHandle);
    glm::mat4 transform = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first

    //change (float)glfwGetTime() to 0.0f to stop the transform
    transform = glm::rotate(transform, (float)glfwGetTime(), glm::vec3(1.0f, 1.0f, 1.0f));

    //get matrix's uniform location and set matrix
    unsigned int transformLoc = glGetUniformLocation(mProgramHandle, "transform");
    glUniformMatrix4fv(transformLoc, 1, GL_FALSE, glm::value_ptr(transform));

    unsigned int viewLoc = glGetUniformLocation(mProgramHandle, "view");
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(cameraView));

    glm::mat4 projection = glm::mat4(1.0f);
    projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    unsigned int projectionLoc = glGetUniformLocation(mProgramHandle, "projection");
    glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));


    
    /*POINT LIGHT*/
    glm::vec3 lightPos = world->pLight.mPoint;
    glm::vec3 pDiffuse = world->pLight.mDiffuse;
    glm::vec3 pSpecular = world->pLight.mSpecular;

    unsigned int lightLoc = glGetUniformLocation(mProgramHandle, "lightPos");
    glUniform3fv(lightLoc, 1, glm::value_ptr(lightPos));

    unsigned int pLightLoc = glGetUniformLocation(mProgramHandle, "pointLight.position");
    glUniform3fv(pLightLoc, 1, glm::value_ptr(lightPos));

    unsigned int pLightDiff = glGetUniformLocation(mProgramHandle, "pointLight.diffuse");
    glUniform3fv(pLightDiff, 1, glm::value_ptr(pDiffuse));

    unsigned int pLightSpec = glGetUniformLocation(mProgramHandle, "pointLight.specular");
    glUniform3fv(pLightSpec, 1, glm::value_ptr(pSpecular));



    /*DIRECTIONAL LIGHT*/
    glm::vec3 direction = world->dLight.mDirection;
    glm::vec3 diffuse = world->dLight.mDiffuse;
    glm::vec3 dSpecular = world->dLight.mSpecular;
    if (!directional) {
        diffuse = { 0, 0, 0 };
        dSpecular = { 0, 0, 0 };
    }

    unsigned int dLightLoc = glGetUniformLocation(mProgramHandle, "dirLight.direction");
    glUniform3fv(dLightLoc, 1, glm::value_ptr(direction));

    unsigned int dLightDiff = glGetUniformLocation(mProgramHandle, "dirLight.diffuse");
    glUniform3fv(dLightDiff, 1, glm::value_ptr(diffuse));

    unsigned int dLightSpec = glGetUniformLocation(mProgramHandle, "pointLight.specular");
    glUniform3fv(dLightSpec, 1, glm::value_ptr(dSpecular));


    //Specular strength of the objet
    unsigned int specularLoc = glGetUniformLocation(mProgramHandle, "specularStr");
    if (!specular) {
        glUniform1f(specularLoc, (GLfloat)0.0f);
    }
    else {
        glUniform1f(specularLoc, (GLfloat)SpecularStrength);
    }

    //Camera position for specular
    glm::vec3 cameraPos = world->cam.mPosition;
    unsigned int camLoc = glGetUniformLocation(mProgramHandle, "cameraPos");
    glUniform3fv(camLoc, 1, glm::value_ptr(cameraPos));

    //move the cube with model matrix
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, cubePos);
    unsigned int modelLoc = glGetUniformLocation(mProgramHandle, "model");
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    glBindVertexArray(mVao);
    glDrawArrays(GL_TRIANGLES, 0, 36);
}

void Square::freeGPUData()
{
    glDeleteVertexArrays(1, &mVao);
    glDeleteBuffers(1, &mVbo);
    glDeleteShader(mFragHandle);
    glDeleteShader(mVertHandle);
    glDeleteProgram(mProgramHandle);
}

// ===---------------POINT LIGHT--------------===

PointLight::PointLight()
{}

PointLight::PointLight(glm::vec3 point, glm::vec3 diffuse, glm::vec3 pSpecular)
{
    mPoint = point;
    mDiffuse = diffuse;
    mSpecular = pSpecular;
}

// ===---------------DIRECTIONAL LIGHT--------------===

DirectionalLight::DirectionalLight()
{}

DirectionalLight::DirectionalLight(glm::vec3 direction, glm::vec3 diffuse, glm::vec3 dSpecular) {
    mDirection = direction;
    mDiffuse = diffuse;
    mSpecular = dSpecular;
}

// ===---------------CAMERA--------------===
Camera::Camera()
{}

Camera::Camera(glm::vec3 position, glm::vec3 target)
{
    mPosition = position;
    mTarget = target;
    calculateRest();
}

//Using the position and target figure out the rest of the vertices for view
void Camera::calculateRest()
{
    mDirection = glm::normalize(mPosition - mTarget);

    mWorldUp = { 0.0f, 1.0f, 0.0f };
    mRight = glm::normalize(glm::cross(mWorldUp, mDirection));
    mUp = glm::cross(mDirection, mRight);
}



// ===---------------INPUT--------------===
void processInput([[maybe_unused]] GLFWwindow* window,
                  [[maybe_unused]] int key,
                  [[maybe_unused]] int scancode,
                  [[maybe_unused]] int action,
                  [[maybe_unused]] int mods) {
    if (key == GLFW_KEY_D && action == GLFW_PRESS) {
        directional = !directional;
    }
    if (key == GLFW_KEY_S && action == GLFW_PRESS) {
        specular = !specular;
    }
}

// ===------------IMPLEMENTATIONS-------------===

Program::Program(int width, int height, std::string title) :
    settings{}, callbacks{}, mWindow{ nullptr }
{
    settings.size.width = width;
    settings.size.height = height;
    settings.title = title;

    if (!glx::initializeGLFW(errorCallback))
    {
        throw OpenGLError("Failed to initialize GLFW with error callback");
    }

    mWindow = glx::createGLFWWindow(settings);
    if (mWindow == nullptr)
    {
        throw OpenGLError("Failed to create GLFW Window");
    }

    createGLContext();
}

void Program::run()
{

    DirectionalLight dLight = world->dLight;

    glfwSetKeyCallback(mWindow, processInput);

    glEnable(GL_DEPTH_TEST);
    while (!glfwWindowShouldClose(mWindow))
    {
        //processInput(mWindow);
        int width;
        int height;
        glfwGetFramebufferSize(mWindow, &width, &height);
        // setup the view to be the window's size
        glViewport(0, 0, width, height);
        // tell OpenGL the what color to clear the screen to
        glClearColor(0, 0, 0, 1);
        // actually clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //render every cube
        for (int i = 0; i < world->scene.size(); i++) {
            world->scene[i].render(world->view, world->cubePositions[i]);
        }

        glfwSwapBuffers(mWindow);
        glfwPollEvents();
    }
}

void Program::freeGPUData()
{
    glx::destroyGLFWWindow(mWindow);
    glx::terminateGLFW();
}

void Program::createGLContext()
{
    glx::bindWindowCallbacks(mWindow, callbacks);
    glfwMakeContextCurrent(mWindow);
    glfwSwapInterval(1);

    if (!glx::createGLContext(mWindow, settings.version))
    {
        throw OpenGLError("Failed to create OpenGL context");
    }

    glx::initializeGLCallback(
        glx::ErrorSource::All, glx::ErrorType::All, glx::ErrorSeverity::All);
}

// ===-----------------DRIVER-----------------===

int main()
{
    try
    {
        // clang-format off
        // clang-format on

        Program prog{ SCR_WIDTH, SCR_HEIGHT, "CSC305 Lab 5" };
      
        std::array<float, 324> vertices2 = {
            //Vertices              //Colors             //normals
            -0.5f, -0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f, -1.0f, //triangle 1 start
             0.5f, -0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f, -1.0f,
             0.5f,  0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f, -1.0f, //triangle 1 end
             0.5f,  0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f, -1.0f, //triangle 2 start
            -0.5f,  0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f, -1.0f,
            -0.5f, -0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f, -1.0f,//triangle 2 end
                                                        
            -0.5f, -0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f,  1.0f, //triangle 3 start
             0.5f, -0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f,  1.0f,
             0.5f,  0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f,  1.0f, //triangle 3 end
             0.5f,  0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f,  1.0f, //triangle 4 start
            -0.5f,  0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f,  1.0f,
            -0.5f, -0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  0.0f,  1.0f, //triangle 4 end
                                                        
            -0.5f,  0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  -1.0f,  0.0f,  0.0f, //triangle 5 start
            -0.5f,  0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  -1.0f,  0.0f,  0.0f,
            -0.5f, -0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  -1.0f,  0.0f,  0.0f, //triangle 5 end
            -0.5f, -0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  -1.0f,  0.0f,  0.0f, //triangle 6 start
            -0.5f, -0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  -1.0f,  0.0f,  0.0f,
            -0.5f,  0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  -1.0f,  0.0f,  0.0f, //triangle 6 end
                                                        
             0.5f,  0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  1.0f,  0.0f,  0.0f, //triangle 7 start
             0.5f,  0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  1.0f,  0.0f,  0.0f,
             0.5f, -0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  1.0f,  0.0f,  0.0f, //triangle 7 end
             0.5f, -0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  1.0f,  0.0f,  0.0f, //triangle 8 start
             0.5f, -0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  1.0f,  0.0f,  0.0f,
             0.5f,  0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  1.0f,  0.0f,  0.0f, //triangle 8 end
                                                        
            -0.5f, -0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f, -1.0f,  0.0f, //triangle 9 start
             0.5f, -0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f, -1.0f,  0.0f,
             0.5f, -0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f, -1.0f,  0.0f, //triangle 9 end
             0.5f, -0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f, -1.0f,  0.0f, //triangle 10 start
            -0.5f, -0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f, -1.0f,  0.0f,
            -0.5f, -0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f, -1.0f,  0.0f, //triangle 10 end
                                                        
            -0.5f,  0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  1.0f,  0.0f, //triangle 11 start
             0.5f,  0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  1.0f,  0.0f,
             0.5f,  0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  1.0f,  0.0f, //triangle 11 end
             0.5f,  0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  1.0f,  0.0f, //triangle 12 start
            -0.5f,  0.5f,  0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  1.0f,  0.0f,
            -0.5f,  0.5f, -0.5f,    1.0f, 0.5f, 0.31f,  0.0f,  1.0f,  0.0f  //triangle 12 end
        };

        Square sqr{};
        Square sqr2{};
        //Declare the point light. Position, diffuse light, specular light
        PointLight point({ 0.0f, 0.0f, 3.0f }, { 1.0f, 1.0f, 1.0f }, {1.0f, 1.0f, 1.0f});
        //Declare directional light. Direction, diffuse light, specular light
        DirectionalLight dir({ 0.0f, 0.0f, 1.0f }, { 1.0f, 1.0f, 1.0f }, { 1.0f, 1.0f, 1.0f });

        //Declare a world, because I really liked it in assignment2

        //Declare where the cubes should be
        glm::vec3 cubePositions[] = {
            glm::vec3(0.0f,  0.0f,  0.0f),
            glm::vec3(-2.0f,  3.0f, 3.0f),
        };

        //Load up both the cubes
        sqr.loadShaders();
        sqr.loadDataToGPU(vertices2);
        sqr.SpecularStrength = 0.0f;
        sqr2.loadShaders();
        sqr2.SpecularStrength = 0.5f;
        sqr2.loadDataToGPU(vertices2);

        //let the world know where the cubes should be
        world->cubePositions.push_back(cubePositions[0]);
        world->cubePositions.push_back(cubePositions[1]);

        //attach both cubes to the scene
        world->scene.push_back(sqr2);
        world->scene.push_back(sqr);

        //attach the point and directional light
        world->pLight = point;
        world->dLight = dir;
        
        //create the camera, set the view and let the world know the camera. (needed for specular)
        Camera camera({ 0, 0, -10 }, { 0, 0, 0 });
        world->view = glm::lookAt(camera.mPosition, camera.mTarget, camera.mUp);
        world->cam = camera;
        prog.run();

        prog.freeGPUData();
        sqr.freeGPUData();
    }
    catch (OpenGLError & err)
    {
        fmt::print("OpenGL Error:\n\t{}\n", err.what());
    }

    return 0;
}
