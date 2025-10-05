from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
#from typing import Protocol, Tuple
# Create and discretize a geometric domain


class Element1D:
    def __init__(self, x_center: float, y_center: float, id: int):
        self.x_center = x_center
        self.y_center = y_center
        self.id = id


class Discretization(ABC):
    @property
    @abstractmethod
    def elements(self) -> list[Element1D]:
        pass

    @property
    @abstractmethod
    def N(self) -> int:
        pass

    @property
    @abstractmethod
    def delta_l(self) -> float:
        pass

    @abstractmethod
    def __getitem__(self, item: int) -> Element1D:
        pass


class Domain(ABC):
    @property
    @abstractmethod
    def geometric_dimensions(self) -> int:
        pass

    @property
    @abstractmethod
    def discretization(self) -> Discretization:
        pass

    @abstractmethod
    def discretize(self):
        pass


class NonUniformElement1D(Element1D):
    def __init__(self, x_center: float, y_center: float, id: int, delta_l):
        super().__init__(x_center, y_center, id)
        self.delta_l = delta_l


class UniformDiscretization(Discretization):
    def __init__(self, delta_l: float, elements: list[Element1D] = None):
        if elements is None:
            self._elements = []
        else:
            self._elements: list[Element1D] = sorted(elements, key=lambda element: element.id)
        self._delta_l = delta_l

    def __add__(self, other):
        if isinstance(other, tuple) and len(other) == 2 and isinstance(other[0], float) and isinstance(other[1], float):
            self.elements.append(Element1D(other[0], other[1], len(self.elements)))
        elif isinstance(other, list):
            for item in other:
                if not isinstance(item, tuple) and len(item) == 2 and isinstance(item[0], float) and isinstance(item[1], float):
                    raise TypeError('List of elements to add must contain only tuples of (x, y) pairs')
                else:
                    self.elements.append(Element1D(item[0], item[1], len(self.elements)))
        else:
            raise TypeError('Not a list of \'Element1D\' objects nor an \'Element1D\' object')
        self.elements = sorted(self.elements, key=lambda element: element.id)

    @property
    def elements(self) -> list[Element1D]:
        return self._elements

    @elements.setter
    def elements(self, elements: list[Element1D]):
        self._elements = elements

    @property
    def delta_l(self) -> float:
        return self._delta_l

    @delta_l.setter
    def delta_l(self, delta_l: float):
        self._delta_l = delta_l

    @property
    def N(self) -> int:
        return len(self.elements)

    def __getitem__(self, index):
        return self.elements[index]

    def __len__(self):
        return len(self.elements)

    @staticmethod
    def color(_id):
        if _id % 2 == 0:
            return 'red'
        else:
            return 'blue'

    def plot_discretization(self, ax: plt.Axes, color=color):
        for element in self.elements:
            x1, y1 = element.x_center - self.delta_l / 2, element.y_center
            x2, y2 = element.x_center + self.delta_l / 2, element.y_center
            ax.plot([x1, x2], [y1, y2], color='red', linestyle='--')
        return ax


class InfiniteMicroStrip(Domain):
    def __init__(self, _width: float, _height: float, _N: int = 10):
        self.w = _width
        self.h = _height
        self.N = _N
        self.delta_l = _width / self.N
        self._discretization = self.discretize()

    @property
    def geometric_dimensions(self):
        return 1

    @property
    def discretization(self) -> Discretization:
        return self._discretization

    @discretization.setter
    def discretization(self, discretization: Discretization):
        self._discretization = discretization

    def discretize(self):
        """
        Use image theorem to reflect the microstrip about the ground plane and discretize both the microstrip and it's reflection

        Parameters
        -------
        self : InfinteMicrostrip(Domain)
            contains all relevant geometric information

        Returns
        -------
        discretization : UniformDiscretization
            a uniform discretization of the domain
        """
        elements = []
        # microstrip
        for i in range(self.N):
            x_coord = i * self.delta_l + self.delta_l / 2  # x_coord = index*delta_l + delta_l / 2
            y_coord = 0
            elements.append(Element1D(x_coord, y_coord, i))

        return UniformDiscretization(delta_l=self.delta_l, elements=elements)


class MicroStrip(Domain):
    def __init__(self, _w_strip: float, _h_strip: float, _w_ground_plane: float, N1_: int = 10, _x_offset_microstrip: float = None):
        self.w_strip = _w_strip
        self.h_strip = _h_strip
        self.w_ground_plane = _w_ground_plane
        self.N1 = N1_

        # set element size to width of microstrip / number of microstrip subdivisions (N1)
        self.delta_l = self.w_strip / self.N1

        # check that w_ground_plane is a multiple of delta_l and determine total number of subdivisions 'N'
        self.N = self.validate_inputs() + self.N1

        # check for x_offset of microstrip
        if _x_offset_microstrip is None:
            # microstrip is centered
            self.x_offset_microstrip = (self.w_ground_plane - self.w_strip) / 2
        else:
            self.x_offset_microstrip = _x_offset_microstrip

        self.discretization = self.discretize()

    def validate_inputs(self):
        if int(self.w_ground_plane % self.delta_l) != 0:
            raise ValueError('Ground plane must be divisible by delta_l = (w_strip / N1)')
        else:
            return int(self.w_ground_plane / self.delta_l)

    @property
    def geometric_dimensions(self):
        return 1

    @property
    def discretization(self) -> Discretization:
        return self._discretization

    @discretization.setter
    def discretization(self, d: Discretization):
        self._discretization = d

    def discretize(self):
        """
        Create a UniformDiscretization of a microstrip transmission line and a finite ground plane

        Parameters
        -------
        self : MicroStrip(Domain)
            The relevant geometric information of the microstrip

        Returns
        -------
        discretization : UniformDiscretization
            a uniform discretization of the microstrip
        """
        elements = []
        for i in range(self.N):
            if i < self.N1:
                # element is part of the microstrip
                y_center = self.h_strip  # y_coord is on the strip of height h
                x_center = i * self.delta_l + self.delta_l / 2 + self.x_offset_microstrip  # x coord is index * delta_l + delta_l / 2 + x_offset of microstrip
                new_element = Element1D(x_center, y_center, i)
                elements.append(new_element)
            else:
                # element is part of ground plane
                y_center = 0  # y_coord is 0
                x_center = (i - self.N1) * self.delta_l + self.delta_l / 2   # x_coord is index of ground plane elements * delta_l + delta_l / 2
                new_element = Element1D(x_center, y_center, i)
                elements.append(new_element)
        return UniformDiscretization(delta_l=self.delta_l, elements=elements)
