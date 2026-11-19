#include <vegafem/volumetricMeshLoader.h>
#include <GL/freeglut.h>

#include <vegafem/corotationalLinearFEM.h>
#include <vegafem/corotationalLinearFEMForceModel.h>
#include <vegafem/generateMassMatrix.h>
#include <vegafem/implicitBackwardEulerSparse.h>

#include <vector>
#include <stdexcept>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>

std::unique_ptr<vegafem::VolumetricMesh> mesh;

double centre[3] = {};
double meshSize = 1.0;
double zoom = 1.0;
double angleX = 20.0;
double angleY = -30.0;

bool dragging = false;
int previousX = 0;
int previousY = 0;

// Declaration order ensures the integrator is destroyed first.
std::unique_ptr<vegafem::SparseMatrix> massMatrix;
std::unique_ptr<vegafem::CorotationalLinearFEM> elasticModel;
std::unique_ptr<vegafem::CorotationalLinearFEMForceModel> forceModel;
std::unique_ptr<vegafem::ImplicitBackwardEulerSparse> integrator;

std::vector<int> fixedDOFs;
std::vector<int> loadedVertices;
std::vector<double> externalForces;

bool loadEnabled = false;
bool paused = false;

constexpr double timestep = 1.0 / 120.0;

// Total force distributed across the free end.
// Demonstration value, in the mesh's consistent force units.
double loadMagnitude = 1000.0;
int loadAxis = 1;

void initialiseSimulation()
{
    const int vertexCount = mesh->getNumVertices();
    const int dofCount = 3 * vertexCount;

    double minimum[3];
    double maximum[3];

    for (int axis = 0; axis < 3; ++axis)
        minimum[axis] = maximum[axis] = mesh->getVertex(0)[axis];

    for (int i = 1; i < vertexCount; ++i)
    {
        const auto& point = mesh->getVertex(i);

        for (int axis = 0; axis < 3; ++axis)
        {
            minimum[axis] = std::min(minimum[axis], point[axis]);
            maximum[axis] = std::max(maximum[axis], point[axis]);
        }
    }

    // Treat the longest bounding-box direction as the beam axis.
    int beamAxis = 0;

    for (int axis = 1; axis < 3; ++axis)
    {
        if (maximum[axis] - minimum[axis] >
            maximum[beamAxis] - minimum[beamAxis])
        {
            beamAxis = axis;
        }
    }

    // Apply a transverse load.
    loadAxis = (beamAxis == 1) ? 2 : 1;

    const double length = maximum[beamAxis] - minimum[beamAxis];
    const double tolerance = length * 1e-6;

    for (int vertex = 0; vertex < vertexCount; ++vertex)
    {
        const double position = mesh->getVertex(vertex)[beamAxis];

        // Fix all three displacement components at this end.
        if (position <= minimum[beamAxis] + tolerance)
        {
            fixedDOFs.push_back(3 * vertex);
            fixedDOFs.push_back(3 * vertex + 1);
            fixedDOFs.push_back(3 * vertex + 2);
        }

        if (position >= maximum[beamAxis] - tolerance)
            loadedVertices.push_back(vertex);
    }

    if (fixedDOFs.empty() || loadedVertices.empty())
        throw std::runtime_error("Could not identify beam ends.");

    vegafem::SparseMatrix* generatedMass = nullptr;

    vegafem::GenerateMassMatrix::computeMassMatrix(
        mesh.get(), &generatedMass, true
    );

    massMatrix.reset(generatedMass);

    elasticModel =
        std::make_unique<vegafem::CorotationalLinearFEM>(mesh.get());

    forceModel =
        std::make_unique<vegafem::CorotationalLinearFEMForceModel>(
            elasticModel.get(), 1
        );

    // Implicit backward Euler with modest Rayleigh damping.
    integrator =
        std::make_unique<vegafem::ImplicitBackwardEulerSparse>(
            dofCount,
            timestep,
            massMatrix.get(),
            forceModel.get(),
            static_cast<int>(fixedDOFs.size()),
            fixedDOFs.data(),
            0.1,     // Mass damping coefficient
            0.001,   // Stiffness damping coefficient
            10,      // Maximum nonlinear iterations
            1e-6     // Convergence tolerance
        );

    integrator->ResetToRest();

    externalForces.assign(dofCount, 0.0);

    std::cout << "Fixed vertices: " << fixedDOFs.size() / 3 << '\n';
    std::cout << "Loaded vertices: " << loadedVertices.size() << '\n';
}

void simulationTimer(int)
{
    if (!paused)
    {
        std::fill(externalForces.begin(), externalForces.end(), 0.0);

        if (loadEnabled)
        {
            const double forcePerVertex =
                -loadMagnitude / loadedVertices.size();

            for (int vertex : loadedVertices)
                externalForces[3 * vertex + loadAxis] = forcePerVertex;
        }

        integrator->SetExternalForces(externalForces.data());

        // Two steps advance approximately 1/60 second.
        for (int step = 0; step < 2; ++step)
        {
            if (integrator->DoTimestep() != 0)
            {
                paused = true;
                std::cerr
                    << "Solver failed. Simulation paused; press R to reset.\n";
                break;
            }
        }

        const double* displacement = integrator->Getq();

        for (int i = 0; i < integrator->GetNumDOFs(); ++i)
        {
            if (!std::isfinite(displacement[i]))
            {
                paused = true;
                std::cerr
                    << "Invalid displacement. Press R to reset.\n";
                break;
            }
        }
    }

    glutPostRedisplay();
    glutTimerFunc(16, simulationTimer, 0);
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glRotated(angleX, 1.0, 0.0, 0.0);
    glRotated(angleY, 0.0, 1.0, 0.0);

    const double scale = 2.0 / meshSize;
    glScaled(scale, scale, scale);
    glTranslated(-centre[0], -centre[1], -centre[2]);

    // Draw all six edges of each tetrahedron.
    const int edges[6][2] = {
        {0, 1}, {0, 2}, {0, 3},
        {1, 2}, {1, 3}, {2, 3}
    };

    glColor3f(0.25f, 0.8f, 1.0f);

    glBegin(GL_LINES);

    for (int element = 0; element < mesh->getNumElements(); ++element)
    {
        for (const auto& edge : edges)
        {
            for (int end = 0; end < 2; ++end)
            {
                const int vertex =
                    mesh->getVertexIndex(element, edge[end]);

                const auto& point = mesh->getVertex(vertex);
                const double* displacement = integrator->Getq();

                glVertex3d(
                    point[0] + displacement[3 * vertex],
                    point[1] + displacement[3 * vertex + 1],
                    point[2] + displacement[3 * vertex + 2]
                );
            }
        }
    }

    glEnd();
    glutSwapBuffers();
}

void reshape(int width, int height)
{
    height = std::max(height, 1);
    width = std::max(width, 1);

    glViewport(0, 0, width, height);

    const double aspect = static_cast<double>(width) / height;
    const double extent = 1.5 * zoom;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (aspect >= 1.0)
        glOrtho(-extent * aspect, extent * aspect,
                -extent, extent, -10.0, 10.0);
    else
        glOrtho(-extent, extent,
                -extent / aspect, extent / aspect, -10.0, 10.0);

    glMatrixMode(GL_MODELVIEW);
}

void updateView()
{
    reshape(glutGet(GLUT_WINDOW_WIDTH),
            glutGet(GLUT_WINDOW_HEIGHT));

    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON)
    {
        dragging = (state == GLUT_DOWN);
        previousX = x;
        previousY = y;
    }
}

void motion(int x, int y)
{
    if (!dragging)
        return;

    angleY += (x - previousX) * 0.5;
    angleX += (y - previousY) * 0.5;

    previousX = x;
    previousY = y;

    glutPostRedisplay();
}

void wheel(int, int direction, int, int)
{
    zoom *= (direction > 0) ? 0.9 : 1.1;
    zoom = std::clamp(zoom, 0.1, 10.0);
    updateView();
}

void keyboard(unsigned char key, int, int)
{
    switch (key)
    {
    case 27:
        glutLeaveMainLoop();
        break;

    case 'l':
    case 'L':
        loadEnabled = !loadEnabled;
        std::cout << "Load "
                  << (loadEnabled ? "ON" : "OFF") << '\n';
        break;

    case 'p':
    case 'P':
        paused = !paused;
        std::cout << (paused ? "Paused\n" : "Running\n");
        break;

    case 'r':
    case 'R':
        loadEnabled = false;
        paused = false;

        integrator->SetExternalForcesToZero();
        integrator->ResetToRest();

        angleX = 20.0;
        angleY = -30.0;
        zoom = 1.0;
        updateView();
        break;
    }
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: UMPS.exe path-to-mesh.veg\n";
        return 1;
    }

    try
    {
        mesh.reset(vegafem::VolumetricMeshLoader::load(argv[1]));

        if (!mesh || mesh->getNumVertices() == 0)
        {
            std::cerr << "Could not load a valid mesh.\n";
            return 1;
        }

        if (mesh->getNumElementVertices() != 4)
        {
            std::cerr << "This viewer requires a tetrahedral mesh.\n";
            return 1;
        }

        double minimum[3];
        double maximum[3];

        for (int axis = 0; axis < 3; ++axis)
            minimum[axis] = maximum[axis] = mesh->getVertex(0)[axis];

        for (int i = 1; i < mesh->getNumVertices(); ++i)
        {
            const auto& point = mesh->getVertex(i);

            for (int axis = 0; axis < 3; ++axis)
            {
                minimum[axis] = std::min(minimum[axis], point[axis]);
                maximum[axis] = std::max(maximum[axis], point[axis]);
            }
        }

        meshSize = 0.0;

        for (int axis = 0; axis < 3; ++axis)
        {
            centre[axis] = (minimum[axis] + maximum[axis]) * 0.5;
            meshSize = std::max(meshSize, maximum[axis] - minimum[axis]);
        }

        if (meshSize <= 0.0)
        {
            std::cerr << "Mesh has no spatial extent.\n";
            return 1;
        }
    }
    catch (...)
    {
        std::cerr << "Error loading mesh.\n";
        return 1;
    }

    try
    {
        initialiseSimulation();
    }
    catch (...)
    {
        std::cerr << "Could not initialise the simulation.\n";
        return 1;
    }

    std::cout << "UMPS mesh viewer\n"
              << "Left-drag: rotate\n"
              << "Mouse wheel: zoom\n"
              << "R: reset view\n"
              << "Escape: close\n";

    // Keep the mesh filename out of GLUT's argument processing.
    int glutArgc = 1;
    glutInit(&glutArgc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1100, 750);
    glutCreateWindow("UMPS - Beam Mesh Viewer");

    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
                  GLUT_ACTION_GLUTMAINLOOP_RETURNS);

    glClearColor(0.07f, 0.09f, 0.13f, 1.0f);
    glEnable(GL_DEPTH_TEST);

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutMouseWheelFunc(wheel);
    glutKeyboardFunc(keyboard);

    glutTimerFunc(16, simulationTimer, 0);

    std::cout << "\nSimulation controls:\n"
              << "L: toggle load\n"
              << "P: pause/resume\n"
              << "R: reset simulation and view\n";

    glutMainLoop();

    return 0;
}