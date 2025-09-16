import domain
import matplotlib.pyplot as plt


if __name__ == '__main__':
    inf_micro_strip = domain.InfiniteMicroStrip(10, 5, 10)
    inf_micro_strip_discretization = inf_micro_strip.discretize()

    micro_strip = domain.MicroStrip(10, 5, 20)
    micro_strip_discretization = micro_strip.discretize()

    def color_microstrip(element_id):
        if micro_strip_discretization[element_id].y_center == 0:
        #if element_id < micro_strip.N1:
            return 'green'
        else:
            return 'blue'

    fig, ax = plt.subplots()

    micro_strip_discretization.plot_discretization(ax, color_microstrip)

    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linestyle="--", alpha=0.6)
    plt.show()
