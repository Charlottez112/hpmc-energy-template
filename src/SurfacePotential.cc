// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "SurfacePotential.h"

namespace hoomd
    {
namespace hpmc
    {

SurfacePotential::SurfacePotential(std::shared_ptr<SystemDefinition> sysdef)
    : ExternalPotential(sysdef), m_params(sysdef->getParticleData()->getNTypes())
    {
    }

LongReal SurfacePotential::particleEnergyImplementation(uint64_t timestep,
                                                        unsigned int tag_i,
                                                        unsigned int type_i,
                                                        const vec3<LongReal>& r_i,
                                                        const quat<LongReal>& q_i,
                                                        LongReal charge_i,
                                                        Trial trial)
    {
    // TODO: implement the functional form of the external potential.
    const auto& param = m_params[type_i];

    // Make the orientation constraint a function of z
    LongReal dist_from_substrate = sqrt((r_i.z - param.m_z_substrate) * (r_i.z - param.m_z_substrate));
    LongReal o_epsilon = exp( - dist_from_substrate * param.m_smoothing_factor) * param.m_orientation_epsilon;
    // Rotate vector by quaternion
    vec3<LongReal> rotated_phi_vec = rotate(q_i, m_phi_ref_vec);
    vec3<LongReal> rotated_theta_vec = rotate(q_i, m_theta_ref_vec);

    // Project phi vector to xy plane
    vec3<LongReal> projected_phi_vec(rotated_phi_vec.x,
                                     rotated_phi_vec.y,
                                     0.0);
    vec3<LongReal> normalized_projected_phi_vec = normalize(projected_phi_vec);

    // Calculate theta for constraint on the flipping motion
    LongReal cos_theta = dot(m_theta_ref_vec, rotated_theta_vec);

    // Calculate phi for constraint on rotation in xy plane
    LongReal cos_phi = dot(m_phi_ref_vec, normalized_projected_phi_vec);
    LongReal cos_4_phi = cos(4 * acos(cos_phi));

    // Constrain the flipping motion
    LongReal theta_arg = (abs(cos_theta) - 1) / param.m_theta_sigma;
    LongReal theta_energy = - exp(-theta_arg * theta_arg) * o_epsilon;

    // Constrain rotation in xy plane
    LongReal phi_arg = (cos_4_phi - 1) / param.m_phi_sigma;
    LongReal phi_energy = - exp(-phi_arg * phi_arg) * o_epsilon;

    // Position energy
    LongReal position_epsilon = param.m_max_position_epsilon * (cos_theta + param.m_y_shift);
    LongReal position_energy = - exp( - dist_from_substrate * param.m_position_sigma) * position_epsilon;

    LongReal energy = position_energy + theta_energy + phi_energy;

    return energy;
    }

void SurfacePotential::setParamsPython(const std::string& particle_type, pybind11::dict params)
    {
    unsigned int particle_type_id = m_sysdef->getParticleData()->getTypeByName(particle_type);
    m_params[particle_type_id] = ParamType(params);
    }

pybind11::dict SurfacePotential::getParamsPython(const std::string& particle_type)
    {
    unsigned int particle_type_id = m_sysdef->getParticleData()->getTypeByName(particle_type);
    return m_params[particle_type_id].asDict();
    }

SurfacePotential::ParamType::ParamType(pybind11::dict params)
    {
    pybind11::dict v = params;

    // TODO: unpack per-type quanties from the Python dictionary to the ParamType struct.

    m_z_substrate = v["z_substrate"].cast<LongReal>();
    m_orientation_epsilon = v["orientation_epsilon"].cast<LongReal>();
    m_theta_sigma = v["theta_sigma"].cast<LongReal>();
    m_phi_sigma = v["phi_sigma"].cast<LongReal>();
    m_smoothing_factor = v["smoothing_factor"].cast<LongReal>();
    m_max_position_epsilon = v["max_position_epsilon"].cast<LongReal>();
    m_y_shift = v["y_shift"].cast<LongReal>();
    m_position_sigma = v["position_sigma"].cast<LongReal>();
    }

pybind11::dict SurfacePotential::ParamType::asDict()
    {
    pybind11::dict pydict;
    pydict["z_substrate"] = m_z_substrate;
    pydict["orientation_epsilon"] = m_orientation_epsilon;
    pydict["theta_sigma"] = m_theta_sigma;
    pydict["phi_sigma"] = m_phi_sigma;
    pydict["smoothing_factor"] = m_smoothing_factor;
    pydict["max_position_epsilon"] = m_max_position_epsilon;
    pydict["y_shift"] = m_y_shift;
    pydict["position_sigma"] = m_position_sigma;
    return pydict;
    }

namespace detail
    {
void export_SurfacePotential(pybind11::module& m)
    {
    pybind11::class_<SurfacePotential, ExternalPotential,
                     std::shared_ptr<SurfacePotential>>(m, "SurfacePotential")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>>())
        .def("setParams", &SurfacePotential::setParamsPython)
        .def("getParams", &SurfacePotential::getParamsPython)
        .def_property("phi_ref_vec",
                      &SurfacePotential::getPhiRefVec,
                      &SurfacePotential::setPhiRefVec)
        .def_property("theta_ref_vec",
                      &SurfacePotential::getThetaRefVec,
                      &SurfacePotential::setThetaRefVec);
    }
    } // end namespace detail
    } // end namespace hpmc
    } // end namespace hoomd
