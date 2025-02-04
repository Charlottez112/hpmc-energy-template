// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#pragma once

#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h"
#include "hoomd/hpmc/ExternalPotential.h"
#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace hpmc
    {

/** Example external potential for use with HPMC simulations.

    This example external potential applies a harmonic trap at the center of the simulation box.

    TODO: Rename the "ExampleExternalPotential" to a class name that represents your potential.
    "ExampleExternal" appears many times in C++, CMakeLists, and Python files. Consider using a
    global search and replace tool.
*/
class SurfacePotential : public ExternalPotential
    {
    public:
    SurfacePotential(std::shared_ptr<SystemDefinition> sysdef);
    virtual ~SurfacePotential() { }

    /// Set type-pair-dependent parameters to the potential.
    void setParamsPython(const std::string& particle_type, pybind11::dict params);

    /// Get type-pair-dependent parameters.
    pybind11::dict getParamsPython(const std::string& particle_type);

    vec3<LongReal> getPhiRefVec() const
    {
        return m_phi_ref_vec;
    }

    void setPhiRefVec(pybind11::tuple phi_vec)
    {
        vec3<LongReal> new_phi_vec(phi_vec[0].cast<LongReal>(),
                                   phi_vec[1].cast<LongReal>(),
                                   phi_vec[2].cast<LongReal>());
        m_phi_ref_vec = new_phi_vec;
    }

    vec3<LongReal> getThetaRefVec() const
    {
        return m_theta_ref_vec;
    }

    void setThetaRefVec(pybind11::tuple theta_vec)
    {
        vec3<LongReal> new_theta_vec(theta_vec[0].cast<LongReal>(),
                                     theta_vec[1].cast<LongReal>(),
                                     theta_vec[2].cast<LongReal>());
        m_theta_ref_vec = new_theta_vec;
    }

    protected:
    LongReal particleEnergyImplementation(uint64_t timestep,
                                          unsigned int tag_i,
                                          unsigned int type_i,
                                          const vec3<LongReal>& r_i,
                                          const quat<LongReal>& q_i,
                                          LongReal charge_i,
                                          Trial trial = Trial::None) override;

    /// per-type parameters
    struct ParamType
        {
        ParamType() { }

        /// Construct a parameter set from a dictionary.
        ParamType(pybind11::dict v);

        /// Convert a parameter set to a dictionary.
        pybind11::dict asDict();

        LongReal m_z_substrate;
        LongReal m_bottom;
        LongReal m_orientation_epsilon;
        LongReal m_theta_sigma;
        LongReal m_phi_sigma;
        LongReal m_smoothing_factor;
        LongReal m_max_position_epsilon;
        LongReal m_scale_diff;
        LongReal m_y_shift;
        LongReal m_position_sigma;
        };

    /// Parameters per type.
    std::vector<ParamType> m_params;
    vec3<LongReal> m_phi_ref_vec;
    vec3<LongReal> m_theta_ref_vec;
    };

namespace detail
    {
//! Export the ExampleUpdater class to python
void export_SurfacePotential(pybind11::module& m);

    } // end namespace detail

    } // end namespace hpmc
    } // end namespace hoomd
