//
//  main.cpp
//  3D Object Drawing
//
//  Created by Nazirul Hasan on 4/9/23.
//

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "shader.h"
#include "sphere.h"
#include "Sphere2.h"
#include "camera.h"
#include "basic_camera.h"
#include "pointLight.h"
#include "cube.h"
#include "stb_image.h"

#include <iostream>

using namespace std;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
//void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);

unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax);

//void bed(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether);
void bed(Cube &cube, Shader& lightingShaderWithTexture, Shader & lightingShader, glm::mat4 alTogether);
void room(Cube &cube, Shader & lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether);
//void rightwall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether);
void table(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether);
void drawPlane(Sphere& sphere, Shader& lightingShader, Shader& lightingShaderWithTexture, glm::mat4 alTogether, Cube& cube, float tx, float ty, float tz);
void Fan(Cube& cube, Shader& lightingShader, glm::mat4 alTogether);
void gate(Cube& cube2, Cube& cube3, Shader& lightingShader, Shader& lightingShaderWithTexture, glm::mat4 alTogether, float tx, float ty, float tz);
void runway(Cube& cube2, Cube& cube3, Shader& lightingShader, Shader& lightingShaderWithTexture, glm::mat4 alTogether, float tx, float ty, float tz);
void room(Cube& cube,Cube &cube2, Cube& cube3,Shader& lightingShader, Shader& lightingShaderWithTexture, glm::mat4 alTogether, float tx, float ty, float tz);

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

// modelling transform
float rotateAngle_X = 0.0;
float rotateAngle_Y = 0.0;
float rotateAngle_Z = 0.0;
float rotateAxis_X = 0.0;
float rotateAxis_Y = 0.0;
float rotateAxis_Z = 1.0;
float translate_X = 0.0;
float translate_Y = 0.0;
float translate_Z = 0.0;
float scale_X = 1.0;
float scale_Y = 1.0;
float scale_Z = 1.0;

// camera
Camera camera(glm::vec3(0.0f, 3.0f, 8.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

float eyeX = 0.0, eyeY = 1.0, eyeZ = 3.0;
float lookAtX = 0.0, lookAtY = 0.0, lookAtZ = 0.0;
glm::vec3 V = glm::vec3(0.0f, 1.0f, 0.0f);
BasicCamera basic_camera(eyeX, eyeY, eyeZ, lookAtX, lookAtY, lookAtZ, V);
float tranz = 0.0;
float  trany = 0.0;



// positions of the point lights
glm::vec3 pointLightPositions[] = {
    glm::vec3(3.50f,  5.50f,  -1.0f), /// Spot
    glm::vec3(0.5f,  1.5f,  0.0f),   /// Point
    glm::vec3(3.0f,  11.5f,  0.5f), /// Sun
    //glm::vec3(-1.5f,  -1.5f,  0.0f)
};
PointLight pointlight1(

    pointLightPositions[0].x, pointLightPositions[0].y, pointLightPositions[0].z,  // position
    //1.0f, 1.0f, 1.0f,     // ambient
    //1.0f, 1.0f, 1.0f,      // diffuse
    //1.0f, 1.0f, 1.0f,        // specular
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    1.0f, 1.0f, 1.0f,        // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    1       // light number
);
PointLight pointlight2(

    pointLightPositions[1].x, pointLightPositions[1].y, pointLightPositions[1].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    1.0f, 1.0f, 1.0f,       // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    2       // light number
);

PointLight pointlight3(

    pointLightPositions[2].x, pointLightPositions[2].y, pointLightPositions[2].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    1.0f, 1.0f, 1.0f,   // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    3       // light number
);

// light settings
bool onOffPointToggle = true;
bool onOffSpotToggle = true;
bool onOffDirectToggle = true;
bool ambientToggle = true;
bool diffuseToggle = true;
bool specularToggle = true;

// timing
float deltaTime = 0.0f;    // time between current frame and last frame
float lastFrame = 0.0f;

void useShaderProgram(Shader& lightingShaderWithTexture)
{
    lightingShaderWithTexture.use();
    pointlight1.setUpPointLight(lightingShaderWithTexture);
    // point light 2
    pointlight2.setUpPointLight(lightingShaderWithTexture);
    // point light 3
    pointlight3.setUpPointLight(lightingShaderWithTexture);
}
Cube *tmp;

Sphere2* spher;


class Curve
{
public:
    vector<float> cntrlPoints;
    vector <float> coordinates;
    vector <float> normals;
    vector <int> indices;
    vector <float> vertices;
    const double pi = 3.14159265389;
    const int nt = 40;
    const int ntheta = 20;
    Curve(vector<float>& tmp)
    {
        this->cntrlPoints = tmp;
        this->fishVAO = hollowBezier(cntrlPoints.data(), ((unsigned int)cntrlPoints.size() / 3) - 1);
        cout << cntrlPoints.size() << endl;
        cout << coordinates.size() << endl;
        cout << normals.size() << endl;
        cout << indices.size() << endl;
        cout << vertices.size() << endl;
    }
    ~Curve()
    {
        glDeleteVertexArrays(1, &fishVAO);
        glDeleteVertexArrays(1, &bezierVAO);
        glDeleteBuffers(1, &bezierVBO);
        glDeleteBuffers(1, &bezierEBO);
    }
    void draw(Shader& lightingShader, glm::mat4 model)
    {
        /// Fish
        lightingShader.use();
        lightingShader.setMat4("model", model);
        lightingShader.setVec3("material.ambient", glm::vec3(1.0f, 0.6f, 0.0f));
        lightingShader.setVec3("material.diffuse", glm::vec3(1.0f, 0.6f, 0.0f));
        lightingShader.setVec3("material.specular", glm::vec3(1.0f, 1.0f, 1.0f));
        lightingShader.setFloat("material.shininess", 32.0f);

        glBindVertexArray(fishVAO);
        glDrawElements(GL_TRIANGLES,                    // primitive type
            (unsigned int)indices.size(),          // # of indices
            GL_UNSIGNED_INT,                 // data type
            (void*)0);                       // offset to indices

        // unbind VAO
        glBindVertexArray(0);
        /// End Fish
    }
private:
    unsigned int fishVAO;
    unsigned int bezierVAO;
    unsigned int bezierVBO;
    unsigned int bezierEBO;


    unsigned int drawControlPoints()
    {
        unsigned int controlPointVAO;
        unsigned int controlPointVBO;

        glGenVertexArrays(1, &controlPointVAO);
        glGenBuffers(1, &controlPointVBO);

        glBindVertexArray(controlPointVAO);

        glBindBuffer(GL_ARRAY_BUFFER, controlPointVBO);
        glBufferData(GL_ARRAY_BUFFER, (unsigned int)cntrlPoints.size() * sizeof(float), cntrlPoints.data(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        return controlPointVAO;
    }

    long long nCr(int n, int r)
    {
        if (r > n / 2)
            r = n - r; // because C(n, r) == C(n, n - r)
        long long ans = 1;
        int i;

        for (i = 1; i <= r; i++)
        {
            ans *= n - r + i;
            ans /= i;
        }

        return ans;
    }
    void BezierCurve(double t, float xy[2], GLfloat ctrlpoints[], int L)
    {
        double y = 0;
        double x = 0;
        t = t > 1.0 ? 1.0 : t;
        for (int i = 0; i < L + 1; i++)
        {
            long long ncr = nCr(L, i);
            double oneMinusTpow = pow(1 - t, double(L - i));
            double tPow = pow(t, double(i));
            double coef = oneMinusTpow * tPow * ncr;
            x += coef * ctrlpoints[i * 3];
            y += coef * ctrlpoints[(i * 3) + 1];

        }
        xy[0] = float(x);
        xy[1] = float(y);
    }
    unsigned int hollowBezier(GLfloat ctrlpoints[], int L)
    {
        int i, j;
        float x, y, z, r;                //current coordinates
        float theta;
        float nx, ny, nz, lengthInv;    // vertex normal


        const float dtheta = 2 * pi / ntheta;        //angular step size

        float t = 0;
        float dt = 1.0 / nt;
        float xy[2];

        for (i = 0; i <= nt; ++i)              //step through y
        {
            BezierCurve(t, xy, ctrlpoints, L);
            r = xy[0];
            y = xy[1];
            theta = 0;
            t += dt;
            lengthInv = 1.0 / r;

            for (j = 0; j <= ntheta; ++j)
            {
                double cosa = cos(theta);
                double sina = sin(theta);
                z = r * cosa;
                x = r * sina;

                coordinates.push_back(x);
                coordinates.push_back(y);
                coordinates.push_back(z);

                // normalized vertex normal (nx, ny, nz)
                // center point of the circle (0,y,0)
                nx = (x - 0) * lengthInv;
                ny = (y - y) * lengthInv;
                nz = (z - 0) * lengthInv;

                normals.push_back(nx);
                normals.push_back(ny);
                normals.push_back(nz);

                theta += dtheta;
            }
        }
        // generate index list of triangles
        // k1--k1+1
        // |  / |
        // | /  |
        // k2--k2+1

        int k1, k2;
        for (int i = 0; i < nt; ++i)
        {
            k1 = i * (ntheta + 1);     // beginning of current stack
            k2 = k1 + ntheta + 1;      // beginning of next stack

            for (int j = 0; j < ntheta; ++j, ++k1, ++k2)
            {
                // k1 => k2 => k1+1
                indices.push_back(k1);
                indices.push_back(k2);
                indices.push_back(k1 + 1);

                // k1+1 => k2 => k2+1
                indices.push_back(k1 + 1);
                indices.push_back(k2);
                indices.push_back(k2 + 1);
            }
        }

        size_t count = coordinates.size();
        for (int i = 0; i < count; i += 3)
        {
            //cout << count << ' ' << i + 2 << endl;
            vertices.push_back(coordinates[i]);
            vertices.push_back(coordinates[i + 1]);
            vertices.push_back(coordinates[i + 2]);

            vertices.push_back(normals[i]);
            vertices.push_back(normals[i + 1]);
            vertices.push_back(normals[i + 2]);
        }

        glGenVertexArrays(1, &bezierVAO);
        glBindVertexArray(bezierVAO);

        // create VBO to copy vertex data to VBO
        glGenBuffers(1, &bezierVBO);
        glBindBuffer(GL_ARRAY_BUFFER, bezierVBO);           // for vertex data
        glBufferData(GL_ARRAY_BUFFER,                   // target
            (unsigned int)vertices.size() * sizeof(float), // data size, # of bytes
            vertices.data(),   // ptr to vertex data
            GL_STATIC_DRAW);                   // usage

        // create EBO to copy index data
        glGenBuffers(1, &bezierEBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bezierEBO);   // for index data
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,           // target
            (unsigned int)indices.size() * sizeof(unsigned int),             // data size, # of bytes
            indices.data(),               // ptr to index data
            GL_STATIC_DRAW);                   // usage

        // activate attrib arrays
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);

        // set attrib arrays with stride and offset
        int stride = 24;     // should be 24 bytes
        glVertexAttribPointer(0, 3, GL_FLOAT, false, stride, (void*)0);
        glVertexAttribPointer(1, 3, GL_FLOAT, false, stride, (void*)(sizeof(float) * 3));

        // unbind VAO, VBO and EBO
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        return bezierVAO;
    }

};

vector<float>Fish = {
    /*
-0.0100, 1.9950, 5.1000,
-0.0550, 1.9800, 5.1000,
-0.0950, 1.9350, 5.1000,
-0.1500, 1.8250, 5.1000,
-0.2250, 1.5900, 5.1000,
-0.2550, 1.3450, 5.1000,
-0.2050, 1.1700, 5.1000,
-0.1400, 1.0050, 5.1000,
-0.0400, 0.8600, 5.1000,
0.0400, 0.7300, 5.1000,
0.1300, 0.6350, 5.1000,
0.2400, 0.5050, 5.1000,*/


-1.4050, 0.4100, 5.1000,
-1.2900, 0.4300, 5.1000,
-1.2100, 0.5000, 5.1000,
-1.1100, 0.5750, 5.1000,
-1.0550, 0.6400, 5.1000,
-1.0200, 0.7200, 5.1000,
-0.9950, 0.7850, 5.1000,
-0.9300, 0.9400, 5.1000,
-0.9150, 0.9850, 5.1000,
-0.8850, 1.0400, 5.1000,
-0.8100, 1.1650, 5.1000,
-0.7750, 1.2300, 5.1000,
-0.7450, 1.2800, 5.1000,
-0.7100, 1.3400, 5.1000,
-0.6750, 1.3950, 5.1000,
-0.5850, 1.4900, 5.1000,
-0.5650, 1.5150, 5.1000,
-0.5350, 1.5350, 5.1000,
-0.4850, 1.5500, 5.1000,
-0.4300, 1.5700, 5.1000,
-0.4100, 1.5800, 5.1000,
-0.3850, 1.5850, 5.1000,
-0.3000, 1.5950, 5.1000,
-0.1950, 1.6050, 5.1000,
-0.1750, 1.6050, 5.1000,
-0.1450, 1.6050, 5.1000,
-0.1200, 1.6050, 5.1000,
-0.0850, 1.6000, 5.1000,
-0.0500, 1.5850, 5.1000,


};
//Cube* tmp, * roomwindow, * roomfloor, * grass, * roomdoor, * walltex, * pondtex, * road1, * divider, * khaja;
Curve* fis;



int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "CSE 4208: Computer Graphics Laboratory", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    //glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile our shader zprogram
    // ------------------------------------
    
    const int parts = 50;
    const float pi = 3.1415926535;
    const float angle = pi * 2.0f / parts;
    float points[200000]{}, radius = 1.0f;

    int ind = 0;
    points[ind++] = 0.0f;
    points[ind++] = 0.0f;
    points[ind++] = 0.0f;
    for (int i = 1; i <= parts; i++) {
        points[ind++] = radius * sin(angle * i);
        points[ind++] = -radius * cos(angle * i);
        points[ind++] = 0.0f;
    }

    for (float r = radius - 0.005f, z = 0.005f; r > 0.0f; r -= 0.005f, z += 0.005f)
    {
        for (int i = 1; i <= parts + 1; i++) {
            points[ind++] = (r + 0.005) * sin(angle * i);
            points[ind++] = -(r + 0.005) * cos(angle * i);
            points[ind++] = z - 0.005f;

            points[ind++] = r * sin(angle * i);
            points[ind++] = -r * cos(angle * i);
            points[ind++] = z;
        }
    }
    for (float r = radius - 0.005f, z = -0.005f; r > 0.0f; r -= 0.005f, z -= 0.005f)
    {
        for (int i = 1; i <= parts + 1; i++) {
            points[ind++] = (r + 0.005) * sin(angle * i);
            points[ind++] = -(r + 0.005) * cos(angle * i);
            points[ind++] = z + 0.005f;

            points[ind++] = r * sin(angle * i);
            points[ind++] = -r * cos(angle * i);
            points[ind++] = z;
        }
    }

    /// Sphere
    unsigned int VBOCL, shpareVAO;
    glGenVertexArrays(1, &shpareVAO);
    glGenBuffers(1, &VBOCL);
    glBindVertexArray(shpareVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBOCL);
    glBufferData(GL_ARRAY_BUFFER, sizeof(points), points, GL_STATIC_DRAW);
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    Shader lightingShaderWithTexture("vertexShaderForPhongShadingWithTexture.vs", "fragmentShaderForPhongShadingWithTexture.fs");
    Shader lightingShader("vertexShaderForPhongShading.vs", "fragmentShaderForPhongShading.fs");
    Shader ourShader("vertexShader.vs", "fragmentShader.fs");

    string diffuseMapPath = "building.png";
    string specularMapPath = "building.png";


    unsigned int diffMap = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    //Cube cube = Cube(diffMap, specMap, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    Cube cube = Cube(diffMap, specMap, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);


    diffuseMapPath = "building.png";
    specularMapPath = "building.png";


    unsigned int diffMap2 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap2 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube2 = Cube(diffMap2, specMap2, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    tmp = &cube2;

    diffuseMapPath = "road.png";
    specularMapPath = "road.png";

    unsigned int diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube3 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    tmp = &cube3;

    diffuseMapPath = "tv.png";
    specularMapPath = "tv.png";


    unsigned int diffMap4 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap4 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube4 = Cube(diffMap4, specMap4, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    tmp = &cube4;

    diffuseMapPath = "sofa.jpg";
    specularMapPath = "sofa.jpg";

    unsigned int diffMap6 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap6 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube6 = Cube(diffMap6, specMap6, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    tmp = &cube6;

    diffuseMapPath = "grass.png";
    specularMapPath = "grass.png";

    unsigned int diffMap5 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap5 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube5 = Cube(diffMap5, specMap5, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    tmp = &cube5;

    diffuseMapPath = "tiles.png";
    specularMapPath = "tiles.png";


    unsigned int diffMap11 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap11 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube11 = Cube(diffMap11, specMap11, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    tmp = &cube11;

    diffuseMapPath = "wall.png";
    specularMapPath = "wall.png";


    unsigned int diffMap12 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap12 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube12 = Cube(diffMap12, specMap12, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    tmp = &cube12;

    diffuseMapPath = "chek.png";
    specularMapPath = "chek.png";


    unsigned int diffMap13 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap13 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube13 = Cube(diffMap13, specMap13, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    tmp = &cube13;

    diffuseMapPath = "red.png";
    specularMapPath = "red.png";


    unsigned int diffMap14 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap14 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube14 = Cube(diffMap14, specMap14, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    tmp = &cube14;

    //Sphere sphere = Sphere();

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    diffuseMapPath = "black.png";
    specularMapPath = "black.png";
    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Sphere2 sp(1.0, 36, 18, glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.5f, 0.5f, 0.5f), 32.0f, diffMap3, specMap3, 0, 1, 0, 1);
    sp.setDefaults();
    sp.setTexture(diffMap3, specMap3);
    spher = &sp;

 

    Sphere ssphere = Sphere();

    Curve fish(Fish);
    fis = &fish;

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.53f, 0.81f, 0.98f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // be sure to activate shader when setting uniforms/drawing objects
        lightingShaderWithTexture.use();
        lightingShaderWithTexture.setVec3("viewPos", camera.Position);

        // pass projection matrix to shader (note that in this case it could change every frame)
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        //glm::mat4 projection = glm::ortho(-2.0f, +2.0f, -1.5f, +1.5f, 0.1f, 100.0f);
        lightingShaderWithTexture.setMat4("projection", projection);

        // camera/view transformation
        glm::mat4 view = camera.GetViewMatrix();
        //glm::mat4 view = basic_camera.createViewMatrix();
        glm::mat4 tmp = glm::translate(glm::mat4(1.0), glm::vec3(0, -trany, tranz));
        lightingShaderWithTexture.setMat4("view", view * tmp);

        lightingShader.use();
        lightingShader.setVec3("viewPos", camera.Position);
        lightingShader.setMat4("projection", projection);
        lightingShader.setMat4("view", view* tmp);

        
        //pointlight1.setUpPointLight(lightingShader);
        //pointlight2.setUpPointLight(lightingShader);
        //pointlight3.setUpPointLight(lightingShader);

        // Modelling Transformation
        glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
        glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;


        //cube.drawCube2(lightingShader, model, 0.345, 0.171, 0.475);
        //lightingShaderWithTexture.use();
        //useShaderProgram(lightingShaderWithTexture);
        //cube2.drawCubeWithTexture(lightingShaderWithTexture, model * glm::translate(glm::mat4(1.0f), glm::vec3( 10, 0, 0) ) );
        //cube4.drawCubeWithTexture(lightingShaderWithTexture, model * glm::translate(glm::mat4(1.0f), glm::vec3(10, 0, 0)));

        //chkin
        glm::mat4 translateMatrix13, rotateXMatrix13, rotateYMatrix13, rotateZMatrix13, scaleMatrix13, model13;
        translateMatrix13 = glm::translate(identityMatrix, glm::vec3(10, .1, 16));
        rotateXMatrix13 = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix13 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix13 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix13 = glm::scale(identityMatrix, glm::vec3(1, 1, 1));
        model13 = translateMatrix13 * rotateXMatrix13 * rotateYMatrix13 * rotateZMatrix13 * scaleMatrix13;
        useShaderProgram(lightingShaderWithTexture);
        cube13.drawCubeWithTexture(lightingShaderWithTexture, model13);

        //redcurpet
        translateMatrix13 = glm::translate(identityMatrix, glm::vec3(20, .2, 13));
        rotateXMatrix13 = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix13 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix13 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix13 = glm::scale(identityMatrix, glm::vec3(80, .1, 14));
        model13 = translateMatrix13 * rotateXMatrix13 * rotateYMatrix13 * rotateZMatrix13 * scaleMatrix13;
        useShaderProgram(lightingShaderWithTexture);
        cube14.drawCubeWithTexture(lightingShaderWithTexture, model13);

        //basement
        identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
        glm::mat4 translateMatrix1, rotateXMatrix1, rotateYMatrix1, rotateZMatrix1, scaleMatrix1, model1;
        translateMatrix1 = glm::translate(identityMatrix, glm::vec3(-100, -.2, -200));
        rotateXMatrix1 = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix1 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix1 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix1 = glm::scale(identityMatrix, glm::vec3(400, .1, 400));
        model1 = translateMatrix1 * rotateXMatrix1 * rotateYMatrix1 * rotateZMatrix1 * scaleMatrix1;
        useShaderProgram(lightingShader);
        cube.drawCube2(lightingShader, model1, .25, 0.25, 0.25);


        //curve object
        identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
        glm::mat4 translateMatrix12, rotateXMatrix12, rotateYMatrix12, rotateZMatrix12, scaleMatrix12, model12;
        translateMatrix12 = glm::translate(identityMatrix, glm::vec3(-40, 1, 5));
        rotateXMatrix12 = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix12 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix12 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix12 = glm::scale(identityMatrix, glm::vec3(5, 5, 5));
        model12 = translateMatrix12 * rotateXMatrix12 * rotateYMatrix12 * rotateZMatrix12 * scaleMatrix12;
        fis->draw(lightingShader, model12);


        //building
        
        for (int i = 15; i < 150; i+=15) {


        translateMatrix1 = glm::translate(identityMatrix, glm::vec3(50, -.2, -80-i));
        rotateXMatrix1 = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix1 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix1 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix1 = glm::scale(identityMatrix, glm::vec3(10, 20, 10));
        model1 = translateMatrix1 * rotateXMatrix1 * rotateYMatrix1 * rotateZMatrix1 * scaleMatrix1;
        useShaderProgram(lightingShaderWithTexture);
        cube2.drawCubeWithTexture(lightingShaderWithTexture, model1 );
        }

        for (int i = 15; i < 150; i += 15) {


            translateMatrix1 = glm::translate(identityMatrix, glm::vec3(60-i, -.2, -180 ));
            rotateXMatrix1 = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
            rotateYMatrix1 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
            rotateZMatrix1 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
            scaleMatrix1 = glm::scale(identityMatrix, glm::vec3(10, 20, 10));
            model1 = translateMatrix1 * rotateXMatrix1 * rotateYMatrix1 * rotateZMatrix1 * scaleMatrix1;
            useShaderProgram(lightingShaderWithTexture);
            cube2.drawCubeWithTexture(lightingShaderWithTexture, model1);
        }



        


        //useShaderProgram(lightingShader);
        drawPlane(ssphere,lightingShader,lightingShaderWithTexture, model, cube, -1.5, 1, -20);
        drawPlane(ssphere, lightingShader, lightingShaderWithTexture, model, cube, -15.5, 1, -20);
        drawPlane(ssphere, lightingShader, lightingShaderWithTexture, model, cube, -29.5, 1, -20);
        drawPlane(ssphere, lightingShader, lightingShaderWithTexture, model, cube, -43.5, 1, -20);




         room(cube, cube4,cube6,lightingShader, lightingShaderWithTexture, model, 55, 0, -30);
         room(cube, cube4, cube6, lightingShader, lightingShaderWithTexture, model, 55, 0, -40);
         room(cube, cube4, cube6, lightingShader, lightingShaderWithTexture, model, 55, 0, -50);
         room(cube, cube4, cube6, lightingShader, lightingShaderWithTexture, model, 55, 0, -50);
         room(cube, cube4, cube6, lightingShader, lightingShaderWithTexture, model, 55, 0, -60);
         runway(cube3, cube5, lightingShader, lightingShaderWithTexture, model, -5, -.1, -10);
         runway(cube3, cube5, lightingShader, lightingShaderWithTexture, model, -19, -.1, -10);
         runway(cube3, cube5, lightingShader, lightingShaderWithTexture, model, -33, -.1, -10);
         runway(cube3, cube5, lightingShader, lightingShaderWithTexture, model, -47, -.1, -10);
         runway(cube3, cube5, lightingShader, lightingShaderWithTexture, model, -61, -.1, -10);
         runway(cube3, cube5, lightingShader, lightingShaderWithTexture, model, -75, -.1, -10);

         gate(cube11, cube12, lightingShader, lightingShaderWithTexture, model,55,0,10);


         //Fan(cube,lightingShader, model);





        

        // also draw the lamp object(s)
        ourShader.use();
        ourShader.setMat4("projection", projection);
        ourShader.setMat4("view", view);

        

        glBindVertexArray(shpareVAO);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }





    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}


void drawPlane(Sphere &sphere,  Shader& lightingShader, Shader& lightingShaderWithTexture,glm::mat4 alTogether, Cube &cube, float tx, float ty, float tz) {
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx, ty+trany, tz-tranz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1, 1, 5));
    model = alTogether * translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;

    sphere.drawSphere(lightingShader, model);

    // right Fan
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx + 0.5, ty+trany, tz - 2-tranz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-30.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(3.5, 0.1, 1));
    model = alTogether * translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 0.25, 0.25, 0.25);


    //wheel
    //identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    //glm::mat4 translateMatrix11, rotateXMatrix11, rotateYMatrix11, rotateZMatrix11, scaleMatrix11, model11;
    //translateMatrix11 = glm::translate(identityMatrix, glm::vec3(tx+0.6, ty-.8, tz+2.5));
    //rotateXMatrix11 = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
    //rotateYMatrix11 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
    //rotateZMatrix11 = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
    //scaleMatrix11 = glm::scale(identityMatrix, glm::vec3(.1, .3, .3));
    //model11 = translateMatrix11 * rotateXMatrix11 * rotateYMatrix11 * rotateZMatrix11 * scaleMatrix11;

    //useShaderProgram(lightingShaderWithTexture);
    ////void drawSphereWithTexture(Shader& lightingShaderWithTexture, glm::mat4 model = glm::mat4(1.0f))
    //spher->drawSphereWithTexture(lightingShaderWithTexture, model11);

    //glm::mat4 translateMatrix11, rotateXMatrix11, rotateYMatrix11, rotateZMatrix11, scaleMatrix11, model11;
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx+0.6, ty-.8+trany, tz+2.5-tranz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(.1, .3, .3));
    model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;

    useShaderProgram(lightingShaderWithTexture);
    ////void drawSphereWithTexture(Shader& lightingShaderWithTexture, glm::mat4 model = glm::mat4(1.0f))
    spher->drawSphereWithTexture(lightingShaderWithTexture, model);



    //wheel-2
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx - 0.6, ty - .8+trany, tz + 2.5-tranz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(.1, .3, .3));
    model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;

    useShaderProgram(lightingShaderWithTexture);
    ////void drawSphereWithTexture(Shader& lightingShaderWithTexture, glm::mat4 model = glm::mat4(1.0f))
    spher->drawSphereWithTexture(lightingShaderWithTexture, model);

    //front-wheel
   // translateMatrix11, rotateXMatrix11, rotateYMatrix11, rotateZMatrix11, scaleMatrix11, model11;
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx + 0.3, ty - .8+trany, tz - 3.5-tranz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(.1, .3, .3));
    model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;

    useShaderProgram(lightingShaderWithTexture);
    ////void drawSphereWithTexture(Shader& lightingShaderWithTexture, glm::mat4 model = glm::mat4(1.0f))
    spher->drawSphereWithTexture(lightingShaderWithTexture, model);

    useShaderProgram(lightingShaderWithTexture);
    //void drawSphereWithTexture(Shader& lightingShaderWithTexture, glm::mat4 model = glm::mat4(1.0f))
    spher->drawSphereWithTexture(lightingShaderWithTexture, model);

    
    //left Fan
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx - 3.5, ty+trany, tz - 0.8-tranz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(30.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(3.5, 0.1, 1));
    model = alTogether * translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 0.25, 0.25, 0.25);

    /*
    //right-back fan
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx + 0.35, ty + 0, tz + 4.5));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-30.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0, 0.1, 0.5));
    model = alTogether * translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 1, 0, 0);

    //left-back fan
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx - 1.18, ty + 0, tz + 5.05));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(30.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1, 0.1, 0.5));
    model = alTogether * translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 1, 0, 0);

    //back-fan
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx + 0, ty + 0.3, tz + 4.5));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(60.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(.1, 1, 0.5));
    model = alTogether * translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 1, 0, 0);*/

    //fan
    // Fan
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx + 0, ty + 0.3+trany, tz +5-tranz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1, 1, 1));
    model = alTogether * translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
    //cube.drawCube2(lightingShader, model1, 0.5, 0.5, 0.5);
    Fan(cube, lightingShader, model);

}

void drawBuilding(Cube cube, Shader lightingShaderWithTexture, glm::mat4 alTogether) {
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0, 0, 0));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1, 1, 1));
    model = alTogether * translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;

    cube.drawCube(lightingShaderWithTexture, model);
}


void room(Cube &cube,Cube &cube2, Cube& cube3, Shader& lightingShader, Shader& lightingShaderWithTexture, glm::mat4 alTogether, float tx, float ty, float tz) {
    
    //room floor
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx, ty, tz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(10, 0.1, 10));
    model = alTogether * translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 0.5, 0.5, 0.5);

    //ceil
    //glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix1, rotateXMatrix1, rotateYMatrix1, rotateZMatrix1, scaleMatrix1, model1;
    translateMatrix1 = glm::translate(identityMatrix, glm::vec3(tx, ty + 10, tz + .2));
    rotateXMatrix1 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix1 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix1 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix1 = glm::scale(identityMatrix, glm::vec3(10, 0.1, 10));
    model1 = alTogether * translateMatrix1 * rotateXMatrix1 * rotateYMatrix1 * rotateZMatrix1 * scaleMatrix1;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model1, 0.5, 0.5, 0.5);

    // Fan
    //glm::mat4 translateMatrix1, rotateXMatrix1, rotateYMatrix1, rotateZMatrix1, scaleMatrix1, model1;
    translateMatrix1 = glm::translate(identityMatrix, glm::vec3(tx + 10/2, ty + 10  - 0.2, tz + .2 + 10/2));
    rotateXMatrix1 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix1 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix1 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix1 = glm::scale(identityMatrix, glm::vec3(1, 1, 1));
    model1 = alTogether * translateMatrix1 * rotateXMatrix1 * rotateYMatrix1 * rotateZMatrix1 * scaleMatrix1;
    useShaderProgram(lightingShader);
    //cube.drawCube2(lightingShader, model1, 0.5, 0.5, 0.5);
    Fan(cube, lightingShader, model1);


    //back wall

    glm::mat4 translateMatrix2, rotateXMatrix2, rotateYMatrix2, rotateZMatrix2, scaleMatrix2, model2;
    translateMatrix2 = glm::translate(identityMatrix, glm::vec3(tx + 10, ty, tz));
    rotateXMatrix2 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix2 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix2 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix2 = glm::scale(identityMatrix, glm::vec3(.2, 10, 10));
    model2 = alTogether * translateMatrix2 * rotateXMatrix2 * rotateYMatrix2 * rotateZMatrix2 * scaleMatrix2;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model2, 0.6784, 0.8471, 0.902);
    //drawCube(cubeVAO, lightingShader, model2, 0.6784, 0.8471, 0.902);

    //right wall
    glm::mat4 translateMatrix3, rotateXMatrix3, rotateYMatrix3, rotateZMatrix3, scaleMatrix3, model3;
    translateMatrix3 = glm::translate(identityMatrix, glm::vec3(tx, ty, tz));
    rotateXMatrix3 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix3 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix3 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix3 = glm::scale(identityMatrix, glm::vec3(10, 10, .2));
    model3 = alTogether * translateMatrix3 * rotateXMatrix3 * rotateYMatrix3 * rotateZMatrix3 * scaleMatrix3;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model3, 0.5569, 0.7294, 0.2588);
    //drawCube(cubeVAO, lightingShader, model3, 0.5569, 0.7294, 0.2588);

    //left wall
    glm::mat4 translateMatrix4, rotateXMatrix4, rotateYMatrix4, rotateZMatrix4, scaleMatrix4, model4;
    translateMatrix4 = glm::translate(identityMatrix, glm::vec3(tx, ty, tz + 10));
    rotateXMatrix4 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix4 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix4 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix4 = glm::scale(identityMatrix, glm::vec3(10, 10, .2));
    model4 = alTogether * translateMatrix4 * rotateXMatrix4 * rotateYMatrix4 * rotateZMatrix4 * scaleMatrix4;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model4, 0.5569, 0.7294, 0.2588);
    //drawCube(cubeVAO, lightingShader, model4, 0.5569, 0.7294, 0.2588);

    //tv
    glm::mat4 translateMatrix5, rotateXMatrix5, rotateYMatrix5, rotateZMatrix5, scaleMatrix5, model5;
    translateMatrix5 = glm::translate(identityMatrix, glm::vec3(tx + 9.5, ty + 2, tz + 2.5));
    rotateXMatrix5 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix5 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix5 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix5 = glm::scale(identityMatrix, glm::vec3(.1, 2, 5));
    model5 = alTogether * translateMatrix5 * rotateXMatrix5 * rotateYMatrix5 * rotateZMatrix5 * scaleMatrix5;
    useShaderProgram(lightingShaderWithTexture);
    cube2.drawCubeWithTexture(lightingShaderWithTexture, model5  );
    

    //seat-left
    glm::mat4 translateMatrix6, rotateXMatrix6, rotateYMatrix6, rotateZMatrix6, scaleMatrix6, model6;
    translateMatrix6 = glm::translate(identityMatrix, glm::vec3(tx, ty + .1, tz + .2));
    rotateXMatrix6 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix6 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix6 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix6 = glm::scale(identityMatrix, glm::vec3(5, 1, 1.2));
    model6 = alTogether * translateMatrix6 * rotateXMatrix6 * rotateYMatrix6 * rotateZMatrix6 * scaleMatrix6;
    useShaderProgram(lightingShaderWithTexture);
    cube3.drawCubeWithTexture(lightingShaderWithTexture, model6);

    //seat-right
    glm::mat4 translateMatrix7, rotateXMatrix7, rotateYMatrix7, rotateZMatrix7, scaleMatrix7, model7;
    translateMatrix7 = glm::translate(identityMatrix, glm::vec3(tx, ty + .1, tz + 8.7));
    rotateXMatrix7 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix7 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix7 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix7 = glm::scale(identityMatrix, glm::vec3(5, 1, 1.2));
    model7 = alTogether * translateMatrix7 * rotateXMatrix7 * rotateYMatrix7 * rotateZMatrix7 * scaleMatrix7;
    useShaderProgram(lightingShaderWithTexture);
    cube3.drawCubeWithTexture(lightingShaderWithTexture, model7);



}

void runway(Cube& cube2, Cube& cube3, Shader& lightingShader, Shader& lightingShaderWithTexture, glm::mat4 alTogether, float tx, float ty, float tz)
{

    //road
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix5, rotateXMatrix5, rotateYMatrix5, rotateZMatrix5, scaleMatrix5, model5;
    translateMatrix5 = glm::translate(identityMatrix, glm::vec3(tx, ty, tz-100));
    rotateXMatrix5 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix5 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix5 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix5 = glm::scale(identityMatrix, glm::vec3(7, .02, 100));
    model5 = alTogether * translateMatrix5 * rotateXMatrix5 * rotateYMatrix5 * rotateZMatrix5 * scaleMatrix5;
    useShaderProgram(lightingShaderWithTexture);
    cube2.drawCubeWithTexture(lightingShaderWithTexture, model5);

    //grass
    //glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix6, rotateXMatrix6, rotateYMatrix6, rotateZMatrix6, scaleMatrix6, model6;
    translateMatrix6 = glm::translate(identityMatrix, glm::vec3(tx + 7, ty , tz -100));
    rotateXMatrix6 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix6 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix6 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix6 = glm::scale(identityMatrix, glm::vec3(7, .02, 100));
    model6 = alTogether * translateMatrix6 * rotateXMatrix6 * rotateYMatrix6 * rotateZMatrix6 * scaleMatrix6;
    useShaderProgram(lightingShaderWithTexture);
    cube3.drawCubeWithTexture(lightingShaderWithTexture, model6);
}


void gate( Cube& cube2, Cube& cube3, Shader& lightingShader, Shader& lightingShaderWithTexture, glm::mat4 alTogether, float tx, float ty, float tz)
{

    //left
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix5, rotateXMatrix5, rotateYMatrix5, rotateZMatrix5, scaleMatrix5, model5;
    translateMatrix5 = glm::translate(identityMatrix, glm::vec3(tx, ty, tz + 20));
    rotateXMatrix5 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix5 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix5 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix5 = glm::scale(identityMatrix, glm::vec3(2, 12, 2));
    model5 = alTogether * translateMatrix5 * rotateXMatrix5 * rotateYMatrix5 * rotateZMatrix5 * scaleMatrix5;
    useShaderProgram(lightingShaderWithTexture);
    cube2.drawCubeWithTexture(lightingShaderWithTexture, model5);

    //right
    //glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix6, rotateXMatrix6, rotateYMatrix6, rotateZMatrix6, scaleMatrix6, model6;
    translateMatrix6 = glm::translate(identityMatrix, glm::vec3(tx, ty, tz ));
    rotateXMatrix6 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix6 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix6 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix6 = glm::scale(identityMatrix, glm::vec3(2, 12, 2));
    model6 = alTogether * translateMatrix6 * rotateXMatrix6 * rotateYMatrix6 * rotateZMatrix6 * scaleMatrix6;
    useShaderProgram(lightingShaderWithTexture);
    cube2.drawCubeWithTexture(lightingShaderWithTexture, model6);

    //upper
    glm::mat4 translateMatrix7, rotateXMatrix7, rotateYMatrix7, rotateZMatrix7, scaleMatrix7, model7;
    translateMatrix7 = glm::translate(identityMatrix, glm::vec3(tx, ty+12, tz));
    rotateXMatrix7 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix7 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix7 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix7 = glm::scale(identityMatrix, glm::vec3(2, 2, 22));
    model7 = alTogether * translateMatrix7 * rotateXMatrix7 * rotateYMatrix7 * rotateZMatrix7 * scaleMatrix7;
    useShaderProgram(lightingShaderWithTexture);
    cube2.drawCubeWithTexture(lightingShaderWithTexture, model7);

    //wall
    glm::mat4 translateMatrix8, rotateXMatrix8, rotateYMatrix8, rotateZMatrix8, scaleMatrix8, model8;
    translateMatrix8 = glm::translate(identityMatrix, glm::vec3(tx, ty, -tz-10));
    rotateXMatrix8 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix8 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix8 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix8 = glm::scale(identityMatrix, glm::vec3(2, 10, 30));
    model8 = alTogether * translateMatrix8 * rotateXMatrix8 * rotateYMatrix8 * rotateZMatrix8 * scaleMatrix8;
    useShaderProgram(lightingShaderWithTexture);
    cube3.drawCubeWithTexture(lightingShaderWithTexture, model8);

    //wall-2
    //glm::mat4 translateMatrix8, rotateXMatrix8, rotateYMatrix8, rotateZMatrix8, scaleMatrix8, model8;
    translateMatrix8 = glm::translate(identityMatrix, glm::vec3(tx, ty, tz + 22));
    rotateXMatrix8 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix8 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix8 = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix8 = glm::scale(identityMatrix, glm::vec3(2, 10, 50));
    model8 = alTogether * translateMatrix8 * rotateXMatrix8 * rotateYMatrix8 * rotateZMatrix8 * scaleMatrix8;
    useShaderProgram(lightingShaderWithTexture);
    cube3.drawCubeWithTexture(lightingShaderWithTexture, model8);
}

float rotateFan = 0;

void Fan(Cube& cube, Shader& lightingShader, glm::mat4 alTogether)
{
    float bladel = 1.5;
    float bladew = 0.2;
    float bladeh = 0.01;

    //glm::mat4 modelForSphere = glm::mat4(1.0f);
    //modelForSphere = glm::translate(alTogether, glm::vec3(1.5f, 1.2f, 0.5f));
    //sphere->setColor(glm::vec3(0.5, 0.2, 0.5));
    //sphere->drawSphere(lightingShader, modelForSphere);

    // Center
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.27, 0.3, 0.27));
    translate = glm::translate(model, glm::vec3(-0.67, 0.0, -0.4));
    //translate = glm::translate(model, glm::vec3(tx, ty, tz));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 0.2, 0.1, 0.5);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0.01, 0.0, 0.0));
    glm::mat4 rotateM = glm::rotate(model, glm::radians(45.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateM * scale * translate;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 1.0f, 0.84f, 0.0f);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0, 0.01, 0.0));
    rotateM = glm::rotate(model, glm::radians(165.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateM * scale * translate;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 1.0f, 0.84f, 0.0f);


    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0.01, 0.01, 0.0));
    rotateM = glm::rotate(model, glm::radians(285.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateM * scale * translate;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 1.0f, 0.84f, 0.0f);

    ////leg1
    //model = glm::mat4(1.0f);
    //translate = glm::mat4(1.0f);
    //translate2 = glm::mat4(1.0f);
    //scale = glm::mat4(1.0f);
    //translate2 = glm::translate(model, glm::vec3(0, 0.02, 0));
    //scale = glm::scale(model, glm::vec3(supportlength, -supporthight, supportwidth));
    //translate = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    //model = alTogether * translate2 * scale * translate;
    //cube.drawCube2(lightingShader, model, 0.804, 0.361, 0.361);
}


// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
    if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS)
    {
        rotateFan += 5.0f;
    }


    if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS)
    {
        tranz += 0.2f;
        lookAtZ += 0.5;
        if (tranz >= 40) {
            trany += 0.1f;
            lookAtY += 0.2;
        }
        

    }
    


    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        camera.ProcessKeyboard(FORWARD, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        camera.ProcessKeyboard(BACKWARD, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        camera.ProcessKeyboard(LEFT, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        camera.ProcessKeyboard(RIGHT, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
        camera.ProcessKeyboard(UP, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        camera.ProcessKeyboard(DOWN, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_UP, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_DOWN, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_LEFT, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_RIGHT, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_LEFT, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_RIGHT, deltaTime+0.2);
    }
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
    {
        if (rotateAxis_X) rotateAngle_X -= 0.1;
        else if (rotateAxis_Y) rotateAngle_Y -= 0.1;
        else rotateAngle_Z -= 0.1;
    }
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) translate_Y += 0.001;
    if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) translate_Y -= 0.001;
    if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) translate_X += 0.001;
    if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) translate_X -= 0.001;
    if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) translate_Z += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) translate_Z -= 0.001;
    ////if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) scale_X += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) scale_X -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS) scale_Y += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS) scale_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) scale_Z += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) scale_Z -= 0.001;

    if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
    {
        rotateAngle_X += 0.1;
        rotateAxis_X = 1.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 0.0;
    }
    if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS)
    {
        rotateAngle_Y += 0.1;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 1.0;
        rotateAxis_Z = 0.0;
    }
    if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
    {
        rotateAngle_Z += 0.1;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 1.0;
    }

    /*
    if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS)
    {
        eyeX += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS)
    {
        eyeX -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS)
    {
        eyeZ += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS)
    {
        eyeZ -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
    {
        eyeY += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
    {
        eyeY -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    */
    if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS)
    {
        pointlight2.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS)
    {
        pointlight2.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS)
    {
        pointlight3.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS)
    {
        pointlight3.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS)
    {
        pointlight1.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
    {
        pointlight1.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnAmbientOn();
        if (pointlight2.isOn())
            pointlight2.turnAmbientOn();
        if (pointlight3.isOn())
            pointlight3.turnAmbientOn();
        //pointlight4.turnDiffuseOn();
        //diffuseToggle = !diffuseToggle;
    //}
    }
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnAmbientOff();
        if (pointlight2.isOn())
            pointlight2.turnAmbientOff();
        if (pointlight3.isOn())
            pointlight3.turnAmbientOff();
        //pointlight4.turnDiffuseOff();
        //diffuseToggle = !diffuseToggle;
    //}
    }

    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnDiffuseOn();
        if (pointlight2.isOn())
            pointlight2.turnDiffuseOn();
        if (pointlight3.isOn())
            pointlight3.turnDiffuseOn();
        //pointlight4.turnAmbientOn();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnDiffuseOff();
        if (pointlight2.isOn())
            pointlight2.turnDiffuseOff();
        if (pointlight3.isOn())
            pointlight3.turnDiffuseOff();
        //diffuseToggle = !diffuseToggle;
        //}
    }


    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOn();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOn();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOn();
        //pointlight4.turnSpecularOn();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        /*cout << "1 " << pointlight1.isOn() << endl;
        cout << pointlight2.isOn() << endl;
        cout << pointlight3.isOn() << endl;*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOff();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOff();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOff();
        //pointlight4.turnSpecularOff();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    // if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    // {
    //     if (onOffPointToggle)
    //     {
    //         pointlight1.turnOff();
    //         
    //         onOffPointToggle = false;
    //     }
    //     else
    //     {
    //         pointlight1.turnOn();
    //       
    //         onOffPointToggle = true;
    //     }
    //    // pointlight3.turnOff();
    //    // pointlight4.turnOff();

    // }
    // 

    // if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    // {
    //     
    //     if (onOffSpotToggle)
    //     {
    //        
    //         pointlight2.turnOff();
    //         onOffSpotToggle = false;
    //     }
    //     else
    //     {
    //         pointlight2.turnOn();
    //         onOffSpotToggle = true;
    //     }
    // }

    // if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
    // {

    //     if (onOffDirectToggle)
    //     {

    //         pointlight3.turnOff();
    //         onOffDirectToggle = false;
    //     }
    //     else
    //     {
    //         pointlight3.turnOn();
    //         onOffDirectToggle = true;
    //     }
    // }
    // 
    // if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
    // {
    //     pointlight1.turnAmbientOn();
    //     pointlight2.turnAmbientOn();
    //    // pointlight3.turnAmbientOn();
    //    // pointlight4.turnAmbientOn();
    // }
    // if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
    // {
    //     pointlight1.turnAmbientOff();
    //     pointlight2.turnAmbientOff();
    //   //  pointlight3.turnAmbientOff();
    //   //  pointlight4.turnAmbientOff();
    // }
    // if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS)
    // {
    //     pointlight1.turnDiffuseOn();
    //     pointlight2.turnDiffuseOn();
    //  //   pointlight3.turnDiffuseOn();
    // //    pointlight4.turnDiffuseOn();
    // }
    // if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS)
    // {
    //     pointlight1.turnDiffuseOff();
    //     pointlight2.turnDiffuseOff();
    ////     pointlight3.turnDiffuseOff();
    // //    pointlight4.turnDiffuseOff();
    // }
    // if (glfwGetKey(window, GLFW_KEY_9) == GLFW_PRESS)
    // {
    //     pointlight1.turnSpecularOn();
    //     pointlight2.turnSpecularOn();
    // //    pointlight3.turnSpecularOn();
    // //    pointlight4.turnSpecularOn();
    // }
    // if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
    // {
    //     pointlight1.turnSpecularOff();
    //     pointlight2.turnSpecularOff();
    ////     pointlight3.turnSpecularOff();
    // //    pointlight4.turnDiffuseOff();
    // }
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}


// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}

unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* data = stbi_load(path, &width, &height, &nrComponents, 0);
    if (data)
    {
        GLenum format;
        if (nrComponents == 1)
            format = GL_RED;
        else if (nrComponents == 3)
            format = GL_RGB;
        else if (nrComponents == 4)
            format = GL_RGBA;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, textureWrappingModeS);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, textureWrappingModeT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, textureFilteringModeMin);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, textureFilteringModeMax);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
    }

    return textureID;
}
