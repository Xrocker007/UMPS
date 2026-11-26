#pragma once

#include <vegafem/forceModel.h>
#include <vegafem/volumetricMesh.h>

#include <array>
#include <memory>
#include <vector>

class PlasticForceModel : public vegafem::ForceModel
{
public:
    PlasticForceModel(vegafem::VolumetricMesh& mesh, double yieldStress);
    void GetInternalForce(const double* u, double* force) override;
    void GetTangentStiffnessMatrixTopology(vegafem::SparseMatrix** matrix) override;
    void GetTangentStiffnessMatrix(const double* u, vegafem::SparseMatrix* matrix) override;
    void GetForceAndMatrix(const double* u, double* force, vegafem::SparseMatrix* matrix) override;
    void ResetToZero() override;

    // Advance history once, after the mechanical step has been accepted.
    bool advancePlasticity(const double* u, double dt);

private:
    struct Element
    {
        std::array<int, 12> dofs;
        std::array<std::array<int, 12>, 12> matrixSlots;
        std::array<vegafem::Vec3d, 4> gradients;
        vegafem::Mat3d inversePlastic;
        double volume;
        double lambda;
        double mu;

        vegafem::Mat3d elasticGradient(const double* u) const;
        vegafem::Mat3d elasticStress(const vegafem::Mat3d& elastic) const;
    };

    std::vector<Element> elements;
    std::vector<vegafem::Mat3d> nextInversePlastic;
    std::unique_ptr<vegafem::SparseMatrix> topology;
    double yieldStress;
};
