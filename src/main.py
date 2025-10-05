import domain
import system_constructor
import matplotlib.pyplot as plt


if __name__ == '__main__':
    inf_micro_strip = domain.InfiniteMicroStrip(10, 5, 100)

    micro_strip = domain.MicroStrip(10, 5, 20)

    inf_constructor = system_constructor.PulseBasisSystem(inf_micro_strip, 10, test_func_type='galerkin')

    print(f'solution coefs: {inf_constructor.rho}')
    fig1, ax1 = plt.subplots()
    inf_constructor.plot_solution(ax1)
    plt.show()

