#include "plasticForceModel.h"

#include <vegafem/volumetricMeshENuMaterial.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using vegafem::Mat3d;
using vegafem::SparseMatrix;
using vegafem::SparseMatrixOutline;
using vegafem::Vec3d;

namespace
{
const Mat3d identity(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
constexpr double plasticFlowRate = 5.0; // per second

double trace(const Mat3d& matrix)
{
    return matrix[0][0] + matrix[1][1] + matrix[2][2];
}
}

PlasticForceModel::PlasticForceModel(vegafem::VolumetricMesh& mesh, double yield)
    : elements(mesh.getNumElements()), nextInversePlastic(elements.size()),
      yieldStress(yield)
{
    r = 3 * mesh.getNumVertices();
    SparseMatrixOutline outline(r);
    for (int i = 0; i < r; ++i)
        outline.AddEntry(i, i);

    for (int index = 0; index < mesh.getNumElements(); ++index)
    {
        auto& element = elements[index];
        auto* material = vegafem::downcastENuMaterial(mesh.getElementMaterial(index));
        if (!material)
            throw std::runtime_error("Plasticity requires isotropic ENU materials.");

        element.lambda = material->getLambda();
        element.mu = material->getMu();
        element.inversePlastic = identity;

        const auto origin = mesh.getVertex(mesh.getVertexIndex(index, 0));
        const auto a = mesh.getVertex(mesh.getVertexIndex(index, 1)) - origin;
        const auto b = mesh.getVertex(mesh.getVertexIndex(index, 2)) - origin;
        const auto c = mesh.getVertex(mesh.getVertexIndex(index, 3)) - origin;
        const Mat3d restEdges = trans(Mat3d(a, b, c));
        element.volume = std::abs(det(restEdges)) / 6.0;
        if (!(element.volume > 0.0) || !std::isfinite(element.volume))
            throw std::runtime_error("Plasticity requires nondegenerate tetrahedra.");

        const Mat3d inverseRest = inv(restEdges);
        element.gradients = {
            -1.0 * (inverseRest[0] + inverseRest[1] + inverseRest[2]),
            inverseRest[0], inverseRest[1], inverseRest[2]
        };

        for (int node = 0; node < 4; ++node)
            for (int axis = 0; axis < 3; ++axis)
                element.dofs[3 * node + axis] = 3 * mesh.getVertexIndex(index, node) + axis;

        for (int row : element.dofs)
            for (int column : element.dofs)
                outline.AddEntry(row, column);
    }

    topology = std::make_unique<SparseMatrix>(&outline);
    for (auto& element : elements)
        for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 12; ++j)
                element.matrixSlots[i][j] = topology->GetInverseIndex(element.dofs[i], element.dofs[j]);
}

Mat3d PlasticForceModel::Element::elasticGradient(const double* u) const
{
    Mat3d total = identity;
    for (int node = 0; node < 4; ++node)
        total += tensorProduct(Vec3d(u + dofs[3 * node]), gradients[node]);

    return total * inversePlastic;
}

Mat3d PlasticForceModel::Element::elasticStress(const Mat3d& elastic) const
{
    const Mat3d strain = 0.5 * (trans(elastic) * elastic - identity);
    return (lambda * trace(strain)) * identity + (2.0 * mu) * strain;
}

void PlasticForceModel::GetTangentStiffnessMatrixTopology(SparseMatrix** matrix)
{
    *matrix = new SparseMatrix(*topology);
}

void PlasticForceModel::GetInternalForce(const double* u, double* force)
{
    GetForceAndMatrix(u, force, nullptr);
}

void PlasticForceModel::GetTangentStiffnessMatrix(const double* u, SparseMatrix* matrix)
{
    GetForceAndMatrix(u, nullptr, matrix);
}

void PlasticForceModel::GetForceAndMatrix(const double* u, double* force, SparseMatrix* matrix)
{
    if (force)
        std::fill_n(force, r, 0.0);
    if (matrix)
        matrix->ResetToZero();

    for (const auto& element : elements)
    {
        const Mat3d elastic = element.elasticGradient(u);
        const Mat3d stress = element.elasticStress(elastic);
        const Mat3d piola = elastic * stress;
        std::array<Vec3d, 4> gradients;
        for (int node = 0; node < 4; ++node)
        {
            gradients[node] = trans(element.inversePlastic) * element.gradients[node];
            if (force)
            {
                const Vec3d nodalForce = element.volume * (piola * gradients[node]);
                for (int axis = 0; axis < 3; ++axis)
                    force[element.dofs[3 * node + axis]] += nodalForce[axis];
            }
        }

        if (!matrix)
            continue;

        // Differentiate the Saint Venant–Kirchhoff force with plastic history held fixed.
        for (int j = 0; j < 12; ++j)
        {
            Vec3d direction(0.0);
            direction[j % 3] = 1.0;
            const Mat3d dElastic = tensorProduct(direction, gradients[j / 3]);
            const Mat3d product = trans(elastic) * dElastic;
            const Mat3d dStress = (element.lambda * trace(product)) * identity
                + element.mu * (product + trans(product));
            const Mat3d dPiola = dElastic * stress + elastic * dStress;

            for (int node = 0; node < 4; ++node)
            {
                const Vec3d derivative = element.volume * (dPiola * gradients[node]);
                for (int axis = 0; axis < 3; ++axis)
                {
                    const int i = 3 * node + axis;
                    matrix->AddEntry(element.dofs[i], element.matrixSlots[i][j], derivative[axis]);
                }
            }
        }
    }
}

bool PlasticForceModel::advancePlasticity(const double* u, double dt)
{
    for (std::size_t i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        nextInversePlastic[i] = element.inversePlastic;
        Mat3d elastic = element.elasticGradient(u);
        const double jacobian = det(elastic);
        if (!(jacobian > 1e-8) || !std::isfinite(jacobian))
            return false;

        const Mat3d cauchy = (elastic * element.elasticStress(elastic) * trans(elastic)) / jacobian;
        const Mat3d deviator = cauchy - (trace(cauchy) / 3.0) * identity;
        const double equivalentStress = std::sqrt(1.5 * (
            dot(deviator[0], deviator[0]) + dot(deviator[1], deviator[1])
            + dot(deviator[2], deviator[2])));
        if (!std::isfinite(equivalentStress))
            return false;
        if (equivalentStress <= yieldStress)
            continue;

        // Apply the inverse of the plastic increment directly because the force evaluation uses Fp inverse.
        const double flow = std::clamp(
            dt * plasticFlowRate * (1.0 - yieldStress / equivalentStress), 0.0, 1.0);
        Mat3d left, right;
        Vec3d singularValues;
        vegafem::SVD(elastic, left, singularValues, right);

        if (!(singularValues[0] > 0 && singularValues[1] > 0 && singularValues[2] > 0))
            return false;

        const double meanLogStretch = (std::log(singularValues[0])
            + std::log(singularValues[1]) + std::log(singularValues[2])) / 3.0;
        Mat3d inverseStretch(0.0);
        for (int axis = 0; axis < 3; ++axis)
            inverseStretch[axis][axis] =
                std::exp(flow * (meanLogStretch - std::log(singularValues[axis])));

        auto& next = nextInversePlastic[i];
        next = right * inverseStretch * trans(right) * element.inversePlastic;
        const double plasticVolume = det(next);
        if (!(plasticVolume > 0) || !std::isfinite(plasticVolume))
            return false;
        next /= std::cbrt(plasticVolume);
    }

    for (std::size_t i = 0; i < elements.size(); ++i)
        elements[i].inversePlastic = nextInversePlastic[i];

    return true;
}

void PlasticForceModel::ResetToZero()
{
    for (auto& element : elements)
        element.inversePlastic = identity;
}
