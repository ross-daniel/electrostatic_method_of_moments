import domain
import system_constructor
import matplotlib.pyplot as plt
import emperical_formulas
import numpy as np
from data_recorder import ErrorTableCollector

collector = ErrorTableCollector()

def record_error(method, wh, N, err):
    collector.add(method, wh, N, err)

h = 1
#w_gnd = 100
#w_strip = 10
#N1 = 20

wh_vals = [0.5, 1, 2]
N1_vals = [5, 10, 20, 50, 100]
#N1_vals = [50]
#wh_vals = [2]

def calc_convergence_rates(errors, delta_ls):
    convergence_rates = []
    for i in range(len(errors)):
        if i == len(errors) - 1:
            break
        else:
            cr = np.log2(errors[i] / errors[i + 1]) / np.log2(delta_ls[i] / delta_ls[i + 1])
            convergence_rates.append(cr)
    return convergence_rates

if __name__ == "__main__":

    #micro_strip_domain = domain.MicroStrip(w_strip, h, w_gnd, N1_=N1)
    #micro_strip_system = system_constructor.MicroStripSystem(micro_strip_domain, test_fn_type=1)

    #fig1, ax1 = plt.subplots()
    #fig2, ax2 = plt.subplots()

    #ax1, ax2 = micro_strip_system.plot_solution([ax1, ax2])

    #ax1.title.set_text('charge vs x-position for microstrip')
    #ax2.title.set_text('charge vs x-position for gnd plane')

    #plt.show()

    for wh in wh_vals:
        errors_pm = []
        errors_galerkin = []
        c_primes_pm = []
        c_primes_galerkin = []
        empirical_c_primes = []
        delta_ls = []
        N_vals = []
        for N1 in N1_vals:
            w = h * wh
            w_gnd = w * 10
            microstrip_discretization = domain.MicroStrip(w, h, w_gnd, N1)
            N = microstrip_discretization.N
            N_vals.append(N)
            delta_ls.append(microstrip_discretization.delta_l)
            pm_system = system_constructor.MicroStripSystem(microstrip_discretization,
                                                            test_fn_type=0)
            galerkin_system = system_constructor.MicroStripSystem(microstrip_discretization, test_fn_type=0)

            #fig1, ax1 = plt.subplots()
            #pm_system.plot_solution([ax1])
            #ax1.title.set_text(f'Surface Charge Density (Point Matching) \n for w_gnd/w_strip = 20, w/h = {wh}, N = {N}')
            #plt.legend()
            #plt.show()

            #fig2, ax2 = plt.subplots()
            #galerkin_system.plot_solution([ax2])
            #ax2.title.set_text(f'Surface Charge Density (Galerkin) \n for w_gnd/w_strip = 20, w/h = {wh}, N = {N}')
            #plt.legend()
            #plt.show()

            empirical_c_prime = emperical_formulas.c_prime(w, h, 1)
            pm_c_prime = pm_system.C_prime
            galerkin_c_prime = galerkin_system.C_prime
            print(f'c_prime: {empirical_c_prime}')
            print(f'pm_c_prime: {pm_c_prime}')
            print(f'galerkin_c_prime: {galerkin_c_prime}')

            errors_pm.append(abs((pm_c_prime - empirical_c_prime) / empirical_c_prime))
            errors_galerkin.append(abs((galerkin_c_prime - empirical_c_prime) / empirical_c_prime))

            c_primes_pm.append(pm_c_prime)
            c_primes_galerkin.append(galerkin_c_prime)
            empirical_c_primes.append(empirical_c_prime)

            #print(f'PM Error for w/h = {wh}, N = {N}: {errors_pm[-1]}')
            record_error('PM', wh, N, errors_pm[-1])
            #print(f'Galerkin Error for w/h = {wh}, N = {N}: {errors_galerkin[-1]}')
            record_error('Galerkin', wh, N, errors_galerkin[-1])
        #plt.plot(N_vals, errors_pm, 'b', label='PM Error')
        #plt.plot(N_vals, errors_galerkin, 'r', label='Galerkin Error')
        #plt.title(f'Error vs N for w/h = {wh}')
        #plt.legend()
        #plt.show()
        #convergence_pm = calc_convergence_rates(errors_pm, delta_ls)
        #convergence_galerkin = calc_convergence_rates(errors_galerkin, delta_ls)
        #plt.plot(N_vals, c_primes_pm, 'blue', label='PM C\'')
        #plt.plot(N_vals, c_primes_galerkin, 'orange', label='Galerkin C\'')
        #plt.plot(N_vals, empirical_c_primes, 'black', label='True C\'')
        #plt.title(f'Convergence for w/h = {wh}')
        #plt.xlabel('N')
        #plt.ylabel('C\'[F/m]')
        #plt.legend()
        #plt.show()
        #plt.plot(N_vals[:-1], convergence_pm, 'blue', label='PM convergence rate')
        #plt.plot(N_vals[:-1], convergence_galerkin, 'orange', label='Galerkin convergence rate')
        #plt.title(f'Convergence for w/h = {wh}')
        #plt.xlabel('N')
        #plt.ylabel('Convergence Rate')
        #plt.legend()
        #plt.show()
        print("% PM table:\n")
        print(collector.to_latex("PM"))
        print("\n\n% Galerkin table:\n")
        print(collector.to_latex("Galerkin"))
        collector.save_csvs()

