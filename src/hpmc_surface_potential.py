# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Template HOOMD-blue component."""

import hoomd
import hoomd.conftest
from hoomd.hpmc_surface_potential import _hpmc_surface_potential


class SurfacePotential(hoomd.hpmc.external.External):
    """Surface potential between CD-RhuA proteins and the mica surface.

    Args:
        z_substrate (float):
        orientation_epsilon (float): 
        theta_sigma (float):
        phi_sigma (float):
        phi_ref_vec (length-3 array):
        theta_ref_vec (length-3 array):
        smoothing_factor (float):
        max_position_epsilon (float):
        y_shift (float):
        position_sigma (float):

    """

    _cpp_class_name = 'SurfacePotential'
    _ext_module = _hpmc_surface_potential

    def __init__(self, phi_ref_vec=(0, 1, 0), theta_ref_vec=(0, 0, 1)):
        param_dict = hoomd.data.parameterdicts.ParameterDict(
            phi_ref_vec=(float, float, float), theta_ref_vec=(float, float, float)
        )
        param_dict['phi_ref_vec'] = phi_ref_vec
        param_dict['theta_ref_vec'] = theta_ref_vec
        self._param_dict = param_dict

        params = hoomd.data.typeparam.TypeParameter(
            'params',
            'particle_types',
            hoomd.data.parameterdicts.TypeParameterDict(z_substrate=float, orientation_epsilon=float,
            theta_sigma=float, phi_sigma=float, smoothing_factor=float, max_position_epsilon=float,
            y_shift=float, position_sigma=float, len_keys=1),
        )
        self._add_typeparam(params)
        self.phi_ref_vec = phi_ref_vec
        self.theta_ref_vec = theta_ref_vec