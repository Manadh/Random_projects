import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.stats import linregress

check = input("\nType in the command line which results you want to display:\n\n\
                a) 1st potential - psi, h and e-foldings.\n\
                b) 2nd potential - psi, h and e-foldings.\n\
                c) 3rd potential - psi, h and e-foldings\n\
                d) Plot of ratio between pressure and energy density\n\
                e) Exercise j), comparison of slopes\n\
                f) Run all\n\n\
                .... i.e input either 'a' or 'b' and so forth:  "
            )


class Cosmos:

    def __init__(self, tau_max, psi_end):

    #arrays and some initial conditions
        self.dtau = 0.01
        self.M = int(tau_max/self.dtau) #number of points in our arrays

        self.psi = np.zeros(self.M)
        self.psi_deriv = np.zeros(self.M)
        self.psi_end = psi_end
        self.ln_a_ai = np.zeros(self.M)
        self.h = np.zeros(self.M)
        self.tau = np.linspace(0,tau_max,self.M)
        self.h[0] = 1
        self.index_inflation_over = 0

    def calculating_arrays(self,psi_i,pot,dpot):

        psi, psi_deriv, h, dtau, ln_a_ai = self.psi, self.psi_deriv, self.h, self.dtau, self.ln_a_ai

        psi[0] = psi_i
        psi_deriv[0] = -(1/3) * dpot(psi_i,psi_i)                               #Condition (11) in the project description
        for j in range(self.M-1):
            pot_j = pot(psi[j],psi_i)                                           #potential in the j'th step
            dpot_j = dpot(psi[j],psi_i)                                          #derivative of the potential in the j'th step

            psi_dderiv = -3*h[j]*psi_deriv[j] - dpot_j                          # equation (9) in the project description
            psi_deriv[j+1] = psi_deriv[j] + psi_dderiv*dtau
            psi[j+1] = psi[j] + psi_deriv[j]*dtau
            h[j+1] = np.sqrt((8 *np.pi/3) * ((1/2)* psi_deriv[j]**2 + pot_j))   #equation (10) in the project description

        if self.psi_end != 1:
            index_inflation_over = np.where(psi <= self.psi_end)[0]             #psi = psi_end  marks the end of inflation
            first_end_index = index_inflation_over[0]
            self.index_inflation_over = first_end_index
        else:
            self.index_inflation_over = -1                                      #arbitrary number set for the third potential, where inflation never ends.

        return psi, h

    def calculating_efoldings(self):
        ln_a_ai = cumtrapz(self.h,self.tau,initial = 0)

        total_efoldings = ln_a_ai[self.index_inflation_over]
        return ln_a_ai, total_efoldings

    def pressure_energydens_ratio(self, pot, psi_i):
        p = (1/2) *self.psi_deriv**2 - pot(self.psi, psi_i)
        rho = (1/2) *self.psi_deriv**2 + pot(self.psi, psi_i)

        ratio = p/rho

        plt.plot(self.tau, ratio)
        plt.xlim(0,200)
        plt.xlabel(r"$\tau$", fontsize = 16)
        plt.ylabel(r"$\dfrac{p_{\phi}}{\rho_{\phi c^2}}$", fontsize = 16, rotation = 0, labelpad = 20)
        plt.axvline(self.tau[self.index_inflation_over], linestyle = "--", color = "r")
        plt.grid()
        plt.legend([r"$\dfrac{p_{\phi}}{\rho_{\phi c^2}}$",r"$\tau_{end}$"], prop = {"size": 18})
        plt.show()





"""
First potential
"""
if check == "a" or check == "f":
    first_pot_object = Cosmos(300, psi_end = 1/np.sqrt(2*np.pi))
    tau1 = first_pot_object.tau

    psi_start1 = 3.10
    pot1 = lambda psi, psi_i: (3 / (8 * np.pi)) * (psi / psi_i)**2
    dpot1 = lambda psi, psi_i: (3 / (4 * np.pi)) * psi / psi_i**2

    psi_array1, h_array1 = first_pot_object.calculating_arrays(psi_start1,pot1,dpot1)
    index_inflation_over = first_pot_object.index_inflation_over
    psi_analytic1 = psi_start1 - (1/(4*np.pi*psi_start1))*tau1

    ln_a_ai, e_foldings1 = first_pot_object.calculating_efoldings()
    print ("e-foldings for potential 1: %.3f" % e_foldings1)

    plt.plot(tau1,psi_array1)
    plt.plot(tau1,psi_analytic1)
    plt.axvline(tau1[0], linestyle = "--", color = "k")
    plt.axvline(tau1[first_pot_object.index_inflation_over], linestyle = "--", color = "k")
    plt.xlabel(r"$\tau$", fontsize = 16)
    plt.ylabel(r"$\psi$", fontsize = 16, rotation = 0)
    plt.grid()
    plt.legend(["Numeric", "Analytic", "Period of inflation"], prop={"size":16})
    plt.show()

    plt.plot(tau1,h_array1)
    plt.xlabel(r"$\tau$", fontsize = 16)
    plt.ylabel("$h = \dfrac{H}{H_i}$", rotation = 0, fontsize = 16, labelpad = 25)
    plt.grid()
    plt.show()



"""
Second potential
"""
if check == "b" or check == "f":

    second_pot_object = Cosmos(3000, psi_end = 1/np.sqrt(np.pi))
    tau2 = second_pot_object.tau

    psi_start2 = 4.41
    pot2 = lambda psi, psi_i: (3/(8*np.pi)) * (psi / psi_i)**4
    dpot2 = lambda psi, psi_i: (3/(2*np.pi)) * (psi**3 / psi_i**4)
    psi_array2, h_array2 = second_pot_object.calculating_arrays(psi_start2,pot2,dpot2)
    psi_analytic2 = psi_start2 * np.exp((-1/(2*np.pi*psi_start2**2))*tau2)
    ln_a_ai2, e_foldings2 = second_pot_object.calculating_efoldings()

    print ("e-foldings for potential 2: %.3f" % e_foldings2)


    plt.plot(tau2,psi_array2)
    plt.plot(tau2,psi_analytic2)
    plt.axvline(tau2[second_pot_object.index_inflation_over], linestyle = "--", color = "k")
    plt.xlabel(r"$\tau$", fontsize = 16)
    plt.ylabel(r"$\psi$", fontsize = 16, rotation = 0)
    plt.grid()
    plt.legend(["Numeric", "Analytic", "End of inflation"], prop={"size":16})
    plt.show()

    plt.plot(tau2,h_array2)
    plt.axvline(tau2[second_pot_object.index_inflation_over], linestyle = "--", color = "k", label = "end of inflation")
    plt.xlabel(r"$\tau$", fontsize = 16)
    plt.ylabel("$h = \dfrac{H}{H_i}$", fontsize = 16, rotation = 0,labelpad = 25)
    plt.grid()
    plt.legend(prop={"size":16})
    plt.show()


"""
Third potential
"""
if check == "c" or check == "f":
    third_pot_object = Cosmos(1000, psi_end = 1)
    tau3 = third_pot_object.tau
    psi_start3 = -5
    lambdaEp = 0.01             #to ensure inflation we set lambdaEp << 4 sqrt(pi) which is what we found as a requirement

    pot3 = lambda psi, psi_i: 3/(8*np.pi) *np.exp(lambdaEp * (psi_i - psi))
    dpot3 = lambda psi, psi_i: -3/(8*np.pi)*lambdaEp * np.exp(lambdaEp * (psi_i - psi))

    psi_array3, h_array3 = third_pot_object.calculating_arrays(psi_start3,pot3,dpot3)
    psi_analytic3 = psi_start3 + (2/lambdaEp)*np.log(1+(lambdaEp)**2/(16*np.pi)*tau3)

    plt.title(r"$\lambda E_P = 0.01$")
    plt.plot(tau3,psi_array3, "-b")
    plt.plot(tau3,psi_analytic3, "--r")
    plt.xlabel(r"$\tau$", fontsize = 16)
    plt.ylabel(r"$\psi$", fontsize = 16, rotation = 0)
    plt.grid()
    plt.legend(["Numeric", "Analytic"], prop={"size":16})
    plt.show()

    plt.plot(tau3,h_array3)
    plt.xlabel(r"$\tau$", fontsize = 16)
    plt.ylabel(r"$h = \dfrac{H}{H_i}$", fontsize = 16, rotation = 0, labelpad = 20)
    plt.grid()
    plt.show()


    #Lets try with different values for lambda
    lambdaEp_list = [0.01,1,4*np.sqrt(np.pi), 8*np.sqrt(np.pi)]

    print ("\nNr of e-foldings for different values of lambdaEp:\n")
    for lambdaEp in lambdaEp_list:
        pot3 = lambda psi, psi_i: 3/(8*np.pi) *np.exp(lambdaEp * (psi_i - psi))
        dpot3 = lambda psi, psi_i: -3/(8*np.pi)*lambdaEp * np.exp(lambdaEp * (psi_i - psi))

        psi_array3, h_array3 = third_pot_object.calculating_arrays(psi_start3,pot3,dpot3)
        psi_analytic3 = psi_start3 + (2/lambdaEp)*np.log(1+(lambdaEp)**2/(16*np.pi)*tau3)
        e_foldings3 = third_pot_object.calculating_efoldings()[1]
        print("lambdaEp = %.3f, e-foldings = %.3f" % (lambdaEp, e_foldings3))

"""
lambdaEp = 0.010, e-foldings = 999.007
lambdaEp = 1.000, e-foldings = 152.956
lambdaEp = 7.090, e-foldings = 7.190
lambdaEp = 14.180, e-foldings = 3.109
"""



if check == "d" or check == "f":

    #exercise h)
    first_pot_object = Cosmos(300, psi_end = 1/np.sqrt(2*np.pi))
    psi_start1 = 3.10
    pot1 = lambda psi, psi_i: (3 / (8 * np.pi)) * (psi / psi_i)**2
    dpot1 = lambda psi, psi_i: (3 / (4 * np.pi)) * psi / psi_i**2
    psi_array1, h_array1 = first_pot_object.calculating_arrays(psi_start1,pot1,dpot1)
    first_pot_object.pressure_energydens_ratio(pot1,psi_start1)


if check == "e" or check == "f":
    #exercise j)
    #As we want to study what happens at later times, I create another object for the first potential

    first_pot_object2 = Cosmos(2000, psi_end = 1/np.sqrt(2*np.pi))
    tau1_2 = first_pot_object2.tau

    psi_start1 = 3.10
    pot1 = lambda psi, psi_i: (3 / (8 * np.pi)) * (psi / psi_i)**2
    dpot1 = lambda psi, psi_i: (3 / (4 * np.pi)) * psi / psi_i**2

    psi_array1_2, h_array1_2 = first_pot_object2.calculating_arrays(psi_start1,pot1,dpot1)
    index_inflation_over = first_pot_object2.index_inflation_over

    ln_a_ai1_2, e_foldings1 = first_pot_object2.calculating_efoldings()

    #normalizing to be able to interpret the plots somewhat
    normalized_array = np.exp(ln_a_ai1_2[index_inflation_over:])/np.exp(ln_a_ai1_2[-1])
    sliced_tau = tau1_2[index_inflation_over:]

    #test functions, also normalized
    norm_tau_two_third = sliced_tau**(2/3)/(tau1_2[-1]**(2/3))
    norm_tau_one_half = sliced_tau**(1/2)/(tau1_2[-1]**(1/2))
    norm_tau_five_sixth = sliced_tau**(5/6)/(tau1_2[-1]**(5/6))
    norm_tau_three_fourth = sliced_tau**(3/4)/(tau1_2[-1]**(3/4))


    index = np.where(tau1_2 >= 1500)                #checking the slope of the curves on tau [1500,2000]
    index_start = index[0][0]
    result = linregress(normalized_array[index_start:], sliced_tau[index_start:])
    result2 = linregress(norm_tau_two_third[index_start:], sliced_tau[index_start:])
    result3 = linregress(norm_tau_one_half[index_start:], sliced_tau[index_start:])
    result4 = linregress(norm_tau_five_sixth[index_start:], sliced_tau[index_start:])
    result5 = linregress(norm_tau_three_fourth[index_start:], sliced_tau[index_start:])

    print ("\nSlopes of numerical solution vs. different powers:")
    print("-------------------------------------------------------------")
    print(
        f"{'Numerical:':<18}{result.slope:> 15.3f}\n",
        f"\n{'p = 2/3:':<18}{result2.slope:> 15.3f}",
        f"\n{'p = 1/2:':<18}{result3.slope:>15.3f}",
        f"\n{'p = 5/6':<18}{result4.slope:> 15.3f}",
        f"\n{'p = 3/4':<18}{result5.slope:> 15.3f}",
    )


    plt.plot(sliced_tau,normalized_array, "--r", label = r"$\dfrac{a(\tau)}{a_i}$")
    plt.plot(sliced_tau,norm_tau_two_third, label = "p = 2/3")
    plt.plot(sliced_tau,norm_tau_one_half, label = "p = 1/2")
    plt.plot(sliced_tau,norm_tau_five_sixth, label = "p = 5/6")
    plt.plot(sliced_tau,norm_tau_three_fourth, label = "p = 3/4")
    plt.xlim(1500,2000)
    plt.ylim(0.75,1.00)
    plt.grid()
    plt.xlabel(r"$\tau$",fontsize = 16)
    plt.ylabel(r"$\tau^p$", rotation = 0,fontsize = 16)
    plt.legend(prop={"size":16})
    plt.show()

"""
Output from j)

Numerical:  2493.124
p = 2/3:    2896.340
p = 1/2:    3793.850
p = 5/6     2358.302
p = 3/4     2597.365
"""
