import numpy as np
from abc import ABC, abstractmethod

from matplotlib import pyplot as plt

import integration
import constants

class SystemConstructor:
    @abstractmethod
    def point_matching_integral(self, i, j):
        pass

    @abstractmethod
    def galerkin_integral(self, i, j):
        pass

    @property
    @abstractmethod
    def integration_parameters(self):
        pass


class PulseBasisSystem(SystemConstructor):
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

    def construct_V_multiple(self, V1, V2):
        # TODO: Implement this method
        # possibly change the domain classes so that the discretization is contained within the
        # 'Domain' subclasses and we have access to N1 here
        return self.construct_V_single(V1 - V2)

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
