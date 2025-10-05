import numpy as np

import domain
import system_constructor
import matplotlib.pyplot as plt
import emperical_formulas

wh_vals = [0.5, 1, 2]
#N_vals = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
N_vals = [5, 10, 20, 30, 40, 50, 100, 150]

h = 1

V1 = 0.1


def calc_convergence_rates(errors, delta_ls):
    convergence_rates = []
    for i in range(len(errors)):
        if i == len(errors) - 1:
            break
        else:
            cr = np.log2(errors[i] / errors[i + 1]) / np.log2(delta_ls[i] / delta_ls[i + 1])
            convergence_rates.append(cr)
    return convergence_rates


if __name__ == '__main__':
    for wh in wh_vals:
        errors_pm = []
        errors_galerkin = []
        c_primes_pm = []
        c_primes_galerkin = []
        empirical_c_primes = []
        delta_ls = []
        for N in N_vals:
            w = h * wh
            inf_microstrip_discretization = domain.InfiniteMicroStrip(w, h, N)
            delta_ls.append(inf_microstrip_discretization.delta_l)
            pm_system = system_constructor.PulseBasisSystem(inf_microstrip_discretization, V1, V2=0, test_func_type='point_matching')
            galerkin_system = system_constructor.PulseBasisSystem(inf_microstrip_discretization, V1, V2=0, test_func_type='galerkin')

            #fig1, ax1 = plt.subplots()
            #pm_system.plot_solution(ax1)
            #ax1.title.set_text(f'Surface Charge Density (Point Matching) for w/h = {wh}, N = {N}')
            #plt.show()

            #fig2, ax2 = plt.subplots()
            #galerkin_system.plot_solution(ax2)
            #ax2.title.set_text(f'Surface Charge Density (Galerkin) for w/h = {wh}, N = {N}')
            #plt.show()


            empirical_c_prime = emperical_formulas.c_prime(w,h,1)
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



            print(f'PM Error for w/h = {wh}, N = {N}: {errors_pm[-1]}')
            print(f'Galerkin Error for w/h = {wh}, N = {N}: {errors_galerkin[-1]}')
        plt.plot(N_vals, errors_pm, 'b', label='PM Error')
        plt.plot(N_vals, errors_galerkin, 'r', label='Galerkin Error')
        plt.title(f'Error vs N for w/h = {wh}')
        plt.legend()
        plt.show()
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

