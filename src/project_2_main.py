import domain
import system_constructor
import matplotlib.pyplot as plt
import emperical_formulas

h = 1
w_gnd = 20
w_strip = 5
N1 = 80

if __name__ == "__main__":

    micro_strip_domain = domain.MicroStrip(w_strip, h, w_gnd, N1_=N1)
    micro_strip_system = system_constructor.MicroStripSystem(micro_strip_domain)

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    ax1, ax2 = micro_strip_system.plot_solution([ax1, ax2])

    ax1.title.set_text('charge vs x-position for microstrip')
    ax2.title.set_text('charge vs x-position for gnd plane')

    plt.show()

