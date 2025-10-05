import numpy as np
from abc import ABC, abstractmethod

from matplotlib import pyplot as plt

import integration
import constants
from domain import MicroStrip, Domain, UniformDiscretization

class SystemConstructor:
    @abstractmethod
    def point_matching_integral(self, i, j):
        pass

    @abstractmethod
    def galerkin_integral(self, i, j):
        pass

    @abstractmethod
    def construct_Z(self):
        pass


class MicroStripSystem(SystemConstructor):
    def __init__(self, domain: MicroStrip, test_fn_type: int = 0, basis_fn_type: int = 0, line_voltages: tuple[float, ...] = None):
        """
        Sets up a system of equations: $[Z][\rho_s] = [V]$ and solves for basis function coefficients of surface charge
        this system uses a finite ground plane

        Parameters
        ----------
        domain: Domain
            A domain that contains the geometric information of the specific problem
        test_fn_type: int
            A number corresponding to the type of test function to use (0=point matching, 1=Galerkin)
        basis_fn_type: int
            A number corresponding to the type of basis function to use (0=pulse functions)
        line_voltages: tuple[float, ...]
            A list of floats corresponding to the voltage of each conductor, line_voltages[0] corresponds to the ground plane
        """
        self.domain = domain
        self.test_fn_type = test_fn_type
        self.basis_fn_type = basis_fn_type
        if line_voltages is None:
            self.line_voltages = [0.5, 1.0]
        else:
            self.line_voltages = line_voltages
        self.num_conductors = len(self.line_voltages)  # currently only 2 conductors are supported but this leaves room for expansion
        # ----- First set up Zij -----
        # Determine integration function based on test function type
        if self.test_fn_type == 0:
            # Point Matching
            self.integral_func = self.point_matching_integral
        elif self.test_fn_type == 1:
            # Galerkin
            self.integral_func = self.galerkin_integral
        else:
            raise ValueError("test_func_type must be 0('point_matching') or 1('galerkin')")
        # $Z_{ij} = \int_{\Omega} \int_{\Omega} W_i * g(x, x') dx'$
        # where W_i is the i'th weighting function and g(x, x') is the greens function for the charge induced at a point x' due to the voltage at point x
        # $ \implies Z_{ij} = - \frac{1}{2*\pi*\eps_0} \int_{\Delta l_i} \int_{Delta l_j} ln(r) dl_j dl_i for a pulse basis system
        self.Z = self.construct_Z()
        # construct the excitation vector 'V'
        # Vi = \int_{\Omega} W_i * V_{12}, V_{12} = V1 if i \in MicroStrip, = V0 if i \in GroundPlane
        self.V = self.construct_V()
        # We want V to be in terms of V1 - V0, not one or the other.
        # Thus, we subtract the first row of Z from every other row, and replace the first row with the condition:
        # \int_{\Omega} \rho_s = 0

        self.adjust_system()
        # Finally, solve the system to find coefficients rho_s
        self.rho = self.solve_system()

        # Calculate C' for error comparisons
        self.C_prime = self.calculate_c_prime()

    def point_matching_integral(self, i, j):
        a = self.domain.discretization.delta_l / 2
        x_diff = self.domain.discretization[i].x_center - self.domain.discretization[j].x_center
        y_diff = self.domain.discretization[i].y_center - self.domain.discretization[j].y_center
        return integration.point_matching_integral(a, x_diff, y_diff)

    def galerkin_integral(self, i, j):
        a = self.domain.discretization.delta_l / 2
        x1 = self.domain.discretization[i].x_center - self.domain.discretization[j].x_center - a
        x2 = self.domain.discretization[i].x_center - self.domain.discretization[j].x_center + a
        y_diff = self.domain.discretization[i].y_center - self.domain.discretization[j].y_center
        return integration.galerkin_integral(a, x1, y_diff, x2, y_diff)

    def construct_Z(self):
        N = self.domain.discretization.N
        Z = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                Z[i, j] = (1 / (2 * np.pi * constants.eps_0)) * self.integral_func(i, j)
        return Z

    def construct_V(self):
        """
        Constructs the V vector for a 2 conductor system
        """
        if self.test_fn_type == 0:
            multiplier = 1
        elif self.test_fn_type == 1:
            multiplier = self.domain.discretization.delta_l
        else:
            raise NotImplementedError('test_fn_type must be 0(point_matching) or 1(galerkin)')

        N = self.domain.discretization.N
        V_gnd = self.line_voltages[0]
        V1 = self.line_voltages[1]
        V = np.zeros((N, 1))
        for i in range(N):
            if i < self.domain.N1:
                V[i] = V1 * multiplier
            else:
                V[i] = V_gnd * multiplier
        return V

    def adjust_system(self):

        Z_1 = self.Z[0, :]  # first row of Z matrix
        V1 = self.line_voltages[1]  # first row of V vector
        N = self.domain.discretization.N

        if self.test_fn_type == 0:
            multiplier = 1
        elif self.test_fn_type == 1:
            multiplier = self.domain.discretization.delta_l
        else:
            raise NotImplementedError('test_fn_type must be 0(point_matching) or 1(galerkin)')

        for i in range(N):
            if i == 0:
                continue
            else:
                self.V[i] -= V1 * multiplier
                for j in range(N):
                    self.Z[i, j] -= Z_1[j]
        # Now set the first row
        for j in range(N):
            self.Z[0, j] = self.domain.discretization.delta_l
        self.V[0] = 0

    def solve_system(self):
        return np.linalg.solve(self.Z, self.V)

    def plot_solution(self, axes: list[plt.Axes]):
        if len(axes) == 1:
            ax1 = axes[0]
            ax2 = axes[0]
        elif len(axes) == 2:
            ax1 = axes[0]
            ax2 = axes[1]
        else:
            raise NotImplementedError('plot_solution only supports 1 or 2 axes')
        # plot charge distribution of MicroStrip on axis 1
        # get an array of the center x vals of the microstrip
        x_vals_strip = [element.x_center for element in self.domain.discretization.elements if element.id < self.domain.N1]
        ax1.plot(x_vals_strip, self.rho[:self.domain.N1], label='strip')
        # plot charge distribution on the ground plane
        x_vals_gnd = [element.x_center for element in self.domain.discretization.elements if element.id >= self.domain.N1]
        ax2.plot(x_vals_gnd, self.rho[self.domain.N1:], label='gnd plane')
        return ax1, ax2

    def calculate_c_prime(self):
        q_prime = 0
        for index, r in enumerate(self.rho):
            if index < self.domain.N1:
                q_prime += r * self.domain.discretization.delta_l
            else:
                break
        return (q_prime / abs((self.line_voltages[0] - self.line_voltages[1]))) * 10 ** -1


class MicrostripInfiniteGndPlaneSystem(SystemConstructor):
    def __init__(self, domain, V1, V2=0, test_func_type='point_matching'):
        self.domain = domain
        self.V0 = V1
        # construct the V vector
        if V2 == 0:
            self.V = self.construct_V_single(V1)
        else:
            self.V = self.construct_V_multiple(V1, V2)

        #print(f'V: {self.V}')
        self.test_func_type = test_func_type
        # determine test func type and save it as a function
        if test_func_type == 'point_matching':
            self.integral_func = self.point_matching_integral
        elif test_func_type == 'galerkin':
            self.integral_func = self.galerkin_integral
        else:
            raise ValueError("test_func_type must be 'point_matching' or 'galerkin'")

        # construct the Z matrix
        self.Z = self.construct_Z()

        #print(f'Z: {self.Z}')

        #Find the solution coefficients
        self.rho = self.solve_system()

        # Calculate C_prime
        self.C_prime = self.calculate_c_prime()

    def point_matching_integral(self, i, j):
        a = self.domain.discretization.delta_l / 2
        x_diff = self.domain.discretization[i].x_center - self.domain.discretization[j].x_center
        return integration.point_matching_integral(a, x_diff, 0) - integration.point_matching_integral(a, x_diff, 2 * self.domain.h)

    def galerkin_integral(self, i, j):
        a = self.domain.discretization.delta_l / 2
        x1 = self.domain.discretization[i].x_center - self.domain.discretization[j].x_center - a
        x2 = self.domain.discretization[i].x_center - self.domain.discretization[j].x_center + a
        return integration.galerkin_integral(a, x1, 0, x2, 0) - integration.galerkin_integral(a, x1, 2 * self.domain.h, x2, 2 * self.domain.h)


    def construct_V_single(self, V):
        return np.zeros((len(self.domain.discretization), 1)) + V

    def construct_Z(self):
        N = len(self.domain.discretization)
        Z = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                Z[i, j] = (1 / (2 * np.pi * constants.eps_0)) * (self.integral_func(i, j))
        return Z

    def solve_system(self):
        if self.test_func_type == 'point_matching':
            multiplier = 1
        else:
            multiplier = self.domain.discretization.delta_l
        #print(f'Vmult: {self.V * multiplier}')
        rho = np.linalg.solve(self.Z, self.V * multiplier)
        return rho

    def calculate_c_prime(self):
        q_prime = 0
        for r in self.rho:
            q_prime += r * self.domain.discretization.delta_l
        return q_prime #/ self.V0

    def plot_solution(self, ax: plt.Axes):
        ax.plot(self.rho)
        return ax
