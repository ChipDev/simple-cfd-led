#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib> // For system()
#include <thread>

class Fluid
{
public:
    enum FieldType
    {
        U_FIELD = 0,
        V_FIELD = 1,
        S_FIELD = 2
        // Add other field types if necessary
    };
    Fluid(float density, int numX, int numY, float h)
        : density(density), numX(numX + 2), numY(numY + 2), h(h)
    {
        numCells = this->numX * this->numY;
        u = std::vector<float>(numCells, 0.0f);
        v = std::vector<float>(numCells, 0.0f);
        newU = std::vector<float>(numCells, 0.0f);
        newV = std::vector<float>(numCells, 0.0f);
        p = std::vector<float>(numCells, 0.0f);
        s = std::vector<float>(numCells, 0.0f);
        m = std::vector<float>(numCells, 1.0f); // Initialized to 1.0
        newM = std::vector<float>(numCells, 0.0f);
    }

    void setupBoundaries()
    {
        int n = numY;
        for (int i = 0; i < numX; i++)
        {
            for (int j = 0; j < numY; j++)
            {
                float sValue = 1.0f; // Assuming fluid by default
                if (i == numX - 1 || j == numY - 1 || j == 0 || i == 0)
                    sValue = 0.0f; // Marking solid boundaries

                s[i * n + j] = sValue;
            }
        }
    }

    void integrate(float dt, float gravity)
    {
        for (int i = 1; i < numX; i++)
        {
            for (int j = 1; j < numY - 1; j++)
            {
                if (s[i * numY + j] != 0.0f && s[i * numY + j - 1] != 0.0f)
                    u[i * numY + j] -= gravity * dt;
            }
        }
    }

    void solveIncompressibility(int numIters, float dt)
    {
        int n = numY;
        float cp = density * h / dt;

        for (int iter = 0; iter < numIters; iter++)
        {
            for (int i = 1; i < numX - 1; i++)
            {
                for (int j = 1; j < numY - 1; j++)
                {
                    if (s[i * n + j] == 0.0f)
                        continue;

                    float sx0 = s[(i - 1) * n + j];
                    float sx1 = s[(i + 1) * n + j];
                    float sy0 = s[i * n + j - 1];
                    float sy1 = s[i * n + j + 1];
                    float s_sum = sx0 + sx1 + sy0 + sy1;
                    if (s_sum == 0.0f)
                        continue;

                    float div = u[(i + 1) * n + j] - u[i * n + j] +
                                v[i * n + j + 1] - v[i * n + j];

                    float p = -div / s_sum;
                    p *= overRelaxation; // Assuming overRelaxation is a member variable
                    this->p[i * n + j] += cp * p;

                    u[i * n + j] -= sx0 * p;
                    u[(i + 1) * n + j] += sx1 * p;
                    v[i * n + j] -= sy0 * p;
                    v[i * n + j + 1] += sy1 * p;
                }
            }
        }
    }

    void extrapolate()
    {
        int n = numY;
        for (int i = 0; i < numX; i++)
        {
            u[i * n] = u[i * n + 1];
            u[i * n + numY - 1] = u[i * n + numY - 2];
        }
        for (int j = 0; j < numY; j++)
        {
            v[j] = v[n + j];
            v[(numX - 1) * n + j] = v[(numX - 2) * n + j];
        }
    }

    float sampleField(float x, float y, FieldType field)
    {
        int n = numY;
        float h1 = 1.0f / h;
        float h2 = 0.5f * h;

        x = std::max(std::min(x, numX * h), h);
        y = std::max(std::min(y, numY * h), h);

        float dx = 0.0f;
        float dy = 0.0f;

        // Pointer to the selected vector
        const std::vector<float> *f;

        switch (field)
        {
        case U_FIELD:
            f = &u;
            dy = h2;
            break;
        case V_FIELD:
            f = &v;
            dx = h2;
            break;
        case S_FIELD:
            f = &m;
            dx = h2;
            dy = h2;
            break;
        }

        int x0 = std::min(static_cast<int>((x - dx) * h1), numX - 1);
        float tx = ((x - dx) - x0 * h) * h1;
        int x1 = std::min(x0 + 1, numX - 1);

        int y0 = std::min(static_cast<int>((y - dy) * h1), numY - 1);
        float ty = ((y - dy) - y0 * h) * h1;
        int y1 = std::min(y0 + 1, numY - 1);

        float sx = 1.0f - tx;
        float sy = 1.0f - ty;

        return sx * sy * (*f)[x0 * n + y0] +
               tx * sy * (*f)[x1 * n + y0] +
               tx * ty * (*f)[x1 * n + y1] +
               sx * ty * (*f)[x0 * n + y1];
    }

    float avgU(int i, int j)
    {
        int n = numY;
        float u_avg = (u[i * n + j - 1] + u[i * n + j] +
                       u[(i + 1) * n + j - 1] + u[(i + 1) * n + j]) *
                      0.25f;
        return u_avg;
    }

    float avgV(int i, int j)
    {
        int n = numY;
        float v_avg = (v[(i - 1) * n + j] + v[i * n + j] +
                       v[(i - 1) * n + j + 1] + v[i * n + j + 1]) *
                      0.25f;
        return v_avg;
    }

    void advectVel(float dt)
    {
        newU = u; // Copying current u values to newU
        newV = v; // Copying current v values to newV

        int n = numY;
        float h = this->h;
        float h2 = 0.5f * h;

        for (int i = 1; i < numX; i++)
        {
            for (int j = 1; j < numY; j++)
            {
                // u component
                if (s[i * n + j] != 0.0f && s[(i - 1) * n + j] != 0.0f && j < numY - 1)
                {
                    float x = i * h;
                    float y = j * h + h2;
                    float u_comp = u[i * n + j];
                    float v_comp = avgV(i, j);
                    x -= dt * u_comp;
                    y -= dt * v_comp;
                    u_comp = sampleField(x, y, U_FIELD);
                    newU[i * n + j] = u_comp;
                }

                // v component
                if (s[i * n + j] != 0.0f && s[i * n + j - 1] != 0.0f && i < numX - 1)
                {
                    float x = i * h + h2;
                    float y = j * h;
                    float u_comp = avgU(i, j);
                    float v_comp = v[i * n + j];
                    x -= dt * u_comp;
                    y -= dt * v_comp;
                    v_comp = sampleField(x, y, V_FIELD);
                    newV[i * n + j] = v_comp;
                }
            }
        }

        u = newU; // Update current u values
        v = newV; // Update current v values
    }

    void advectSmoke(float dt)
    {
        newM = m; // Copying current m values to newM

        int n = numY;
        float h = this->h;
        float h2 = 0.5f * h;

        for (int i = 1; i < numX - 1; i++)
        {
            for (int j = 1; j < numY - 1; j++)
            {
                if (s[i * n + j] != 0.0f)
                {
                    float u_comp = (u[i * n + j] + u[(i + 1) * n + j]) * 0.5f;
                    float v_comp = (v[i * n + j] + v[i * n + j + 1]) * 0.5f;
                    float x = i * h + h2 - dt * u_comp;
                    float y = j * h + h2 - dt * v_comp;

                    newM[i * n + j] = sampleField(x, y, S_FIELD);
                }
            }
        }

        m = newM; // Update current m values
    }

    void simulate(float dt, float gravity, int numIters)
    {
        // Integrate external forces like gravity
        integrate(dt, gravity);

        // Reset pressure values to zero
        std::fill(p.begin(), p.end(), 0.0f);

        // Solve the incompressibility condition
        solveIncompressibility(numIters, dt);

        // Extrapolate velocities at the boundaries
        extrapolate();

        // Advect velocity fields
        advectVel(dt);

        // Advect scalar field (like smoke)
        advectSmoke(dt);
    }

    float density;
    int numX, numY, numCells;
    float h;
    std::vector<float> u, v, newU, newV, p, s, m, newM;
    // ... other member variables ...

    float overRelaxation = 1.9f;
};

struct Scene
{
    float gravity;
    float dt;
    int numIters;
    float overRelaxation;
    bool showSmoke;
    Fluid fluid; // Fluid object

    Scene()
        : gravity(-0.0f), dt(1.0f / 60.0f), numIters(2),
          overRelaxation(1.9f), showSmoke(true),
          // Initialize Fluid directly here
          fluid(500.0f, 10, 10, 0.1f) //doenst change
    { // Example parameters for Fluid
      // Additional setup can be done here if needed
    }

    void setupScene(int sceneNr);
};

void Scene::setupScene(int sceneNr)
{

    // Existing initialization code
    float density = 1.0f; // or 1.0f if density is uniform
    float h = 1.0f / 100;     // Assuming grid size of 50 for example
    fluid = Fluid(density, 100, 100, h);
    fluid.setupBoundaries();

    // Add a disturbance in the middle of the grid
    int midX = fluid.numX / 2;
    int midY = fluid.numY / 2;

    // Disturbance in the scalar field (e.g., smoke density)
    fluid.m[midX * fluid.numY + midY] = 1.0f; // Example value

    // Initial velocity (e.g., upward flow)
    for (int i = 1; i < fluid.numX - 1; ++i)
    {
        for (int j = 1; j < fluid.numY - 1; ++j)
        {
            fluid.u[i * fluid.numY + j] = 0.0f; // Horizontal component
            fluid.v[i * fluid.numY + j] = 0.0f; // Vertical component, example value
        }
    }
}

void setColor(float value)
{
    // Clamp the value to the range [0, 5]
    if (value < 0.0f)
        value = 0.0f;
    if (value > 5.0f)
        value = 5.0f;

    // Map the value to the 256-color palette
    int colorIndex;
    if (value < 0.5f)
    {
        // Map the range [0, 0.5] to a smaller part of the palette
        colorIndex = 16 + static_cast<int>((value / 0.5f) * 20); // Example: using first 20 colors
    }
    else if (value <= 1.5f)
    {
        // Emphasize the range [0.5, 1.5] over a larger part of the palette
        colorIndex = 36 + static_cast<int>(((value - 0.5f) / 1.0f) * 150); // Example: using next 150 colors
    }
    else
    {
        // Map the rest of the range [1.5, 5] to the remaining part of the palette
        colorIndex = 186 + static_cast<int>(((value - 1.5f) / 3.5f) * 45); // Example: using last 45 colors
    }

    // Set the color using ANSI escape code
    std::cout << "\033[48;5;" << colorIndex << "m";
}

// Function to generate ANSI escape code for RGB background color
std::string getAnsiBackgroundCode(float r, float g, float b)
{
    int ri = static_cast<int>(r * 255);
    int gi = static_cast<int>(g * 255);
    int bi = static_cast<int>(b * 255);
    return "\033[48;2;" + std::to_string(ri) + ";" + std::to_string(gi) + ";" + std::to_string(bi) + "m";
}

std::string getSciColor(float val, float minVal, float maxVal, float sv)
{
    val = std::min(std::max(val, minVal), maxVal - 0.0001f);
    float d = maxVal - minVal;
    val = (d == 0.0f) ? 0.5f : (val - minVal) / d;
    float m = 0.25f;
    int num = static_cast<int>(val / m);
    float s = (val - num * m) / m;

    float r, g, b;
    switch (num)
    {
    case 0:
        r = 0.0f;
        g = s;
        b = 1.0f;
        break;
    case 1:
        r = 0.0f;
        g = 1.0f;
        b = 1.0f - s;
        break;
    case 2:
        r = s;
        g = 1.0f;
        b = 0.0f;
        break;
    case 3:
        r = 1.0f;
        g = 1.0f - s;
        b = 0.0f;
        break;
    default:
        r = g = b = 0.0f;
        break; // Default case to handle unexpected values
    }

    return getAnsiBackgroundCode(r,g,b);
}

int main()
{
    // ... setup your fluid simulation ...
    Scene scene;
    scene.setupScene(0);

    // Number of simulation steps
    int steps = 5000;
    for (int step = 0; step < steps; ++step)
    {

        int midX = scene.fluid.numX / 2;
        int midY = scene.fluid.numY / 2;
        for (int i = 15; i < 19; ++i)
        {

            if(step < 300) {
                scene.fluid.m[midX * scene.fluid.numY] = 2.0f; // Example value
                scene.fluid.m[(midX + 1) * scene.fluid.numY] = 2.0f; 
                scene.fluid.m[(midX + 2) * scene.fluid.numY] = 2.0f; 
                scene.fluid.m[(midX - 1) * scene.fluid.numY] = 2.0f; 
            }

            scene.fluid.v[midX * scene.fluid.numY + 1] = 1.0f; // Example value
            scene.fluid.v[(midX + 1) * scene.fluid.numY + 1] = 2.0f; 

            scene.fluid.s[(midX-1) * scene.fluid.numY + i] = 0.0f;
            scene.fluid.s[(midX) * scene.fluid.numY + i] = 0.0f;
            scene.fluid.s[(midX+1) * scene.fluid.numY + i] = 0.0f;
            scene.fluid.s[(midX+2) * scene.fluid.numY + i] = 0.0f;
        }

        // ... update the fluid simulation ...
        std::cout << "\033[H";
        for (int i = 0; i < scene.fluid.numX; ++i)
        {
            for (int j = 0; j < scene.fluid.numY; ++j)
            {

                float value = scene.fluid.m[i * scene.fluid.numY + j] * scene.fluid.s[i * scene.fluid.numY + j];
                std::cout << getSciColor(value, 1.0f, 2.0f, scene.fluid.s[i * scene.fluid.numY + j] * scene.fluid.s[i * scene.fluid.numY + j]); // Set color based on value
                bool obst = (scene.fluid.s[i * scene.fluid.numY + j] * scene.fluid.s[i * scene.fluid.numY + j] == 0);
                
                // Print the value with 3 decimal places
                std::cout << std::fixed << std::setprecision(1) << std::setw(1) << (obst?"##":"  ");
            }
            std::cout << "\033[0m\n"; // Reset color and move to new line
        }
        std::cout << "----------------------" << std::endl;
        scene.fluid.simulate(scene.dt, scene.gravity, scene.numIters);
        //std::this_thread::sleep_for(std::chrono::milliseconds(10)); // Adjust delay as needed
    }

    return 0;
}

// int main()
// {
//     Scene scene;
//     scene.setupScene(0);

//     // Number of simulation steps
//     int steps = 5;
//     for (int step = 0; step < steps; ++step)
//     {

//         // Print the current state of the grid (example: printing smoke density)
//         if (scene.showSmoke)
//         {
//             for (int i = 1; i < scene.fluid.numX - 1; ++i)
//             {
//                 for (int j = 1; j < scene.fluid.numY - 1; ++j)
//                 {
//                     std::cout << std::fixed << std::setprecision(2) << scene.fluid.m[i * scene.fluid.numY + j] << " ";
//                 }
//                 std::cout << std::endl;
//             }
//             std::cout << "----------------------" << std::endl;
//         }

//          // Update the fluid simulation
//         scene.fluid.simulate(scene.dt, scene.gravity, scene.numIters);
//     }

//     return 0;
// }