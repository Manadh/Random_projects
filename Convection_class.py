#rho visulaliser

import FVis3
import numpy as np


class The2Dconvection:


    def __init__(self):


        #Pertubation check
        check = input("Pertubation on/off? (y/n) ")

        if check == "y":

            self.Pertubation = True

        elif check == "n":

            self.Pertubation = False

        else:

            print ("Please enter either 'y' or 'n'")

            exit()


        #defining constants
        self.my = 0.61
        self.mu = 1.66E-27           # [kg]
        self.kB = 1.38E-23          # [m^2 kg s^-2 K^-1]
        self.gamma = 5/3            # adiabatic index
        self.P0 = 1.8E8             # [Pa] [N/m^2]
        self.nabla = 0.4101         #Slighly higher than 2/5 to enable convection
        self.G = 6.67E-11           # [m^3 kg^-1 s^-2]
        self.M0 = 1.989E30          # [kg]
        self.R0 = 6.96E8            # [m]
        self.T0 = 5778              #[K]
        self.g = self.G * self.M0 /(self.R0**2)    #we asssume g = constant, as we move inwards by a small amount (~1.7% R_0)

        #Nr of points in our arrays
        self.nx = 300
        self.ny = 100

        #Start and end values for our enviroment [in meters]
        self.x0 = 0
        self.y0 = 0
        self.y_end = 4E6
        self.x_end = 12E6

        #Step length of grid points
        self.dx = self.x_end/(self.nx - 1)
        self.dy = self.y_end/(self.ny - 1)

        #grid arrays
        self.x = np.linspace(self.x0,self.x_end,self.nx)
        self.y = np.linspace(self.y0, self.y_end,self.ny)

        ##matrices
        nx = self.nx
        ny = self.ny

        self.P = np.zeros((ny,nx))
        self.T = np.zeros((ny,nx))

        self.rho = np.zeros((ny,nx))
        self.e = np.zeros((ny,nx))
        self.u = np.zeros((ny,nx))
        self.w = np.zeros((ny,nx))

        self.rho_u = self.rho * self.u
        self.rho_w = self.rho * self.w


        #Temporal derivatives
        self.drhodt = np.zeros((ny,nx))
        self.dedt = np.zeros((ny,nx))
        self.dudt = np.zeros((ny,nx))
        self.dwdt = np.zeros((ny,nx))
        self.dPdt = np.zeros((ny,nx))
        self.drho_udt = np.zeros((ny,nx))
        self.drho_wdt = np.zeros((ny,nx))


        #Pertubation matrix
        self.Pert = np.zeros((ny,nx))
        self.Pert2 = np.zeros((ny,nx))


    def initialise(self):

        """
        initialise temperature, pressure, density and internal energy
        """

        my, mu, kB, gamma, nabla, g = self.my, self.mu, self.kB, self.gamma, self.nabla, self.g


        #Initialise the temperature array

        #Here we set it up it so that the highest value for T happens at y = 0, thus at y = 0, we are at the bottom of
        #our defined space, i.e closer to the core of the star. Since
        for j in range(self.ny):
            self.T[j,:] = self.T0 - (nabla*my*mu*g*(self.y[j]-self.y_end))/kB


        # update P early to get the circular arrow plot

        self.P = self.P0*(self.T/self.T0)**(1/nabla)

        #Here the pertubation is activated if it was enabled by user
        if self.Pertubation == True:

            """
            one central bomb
            """
            sigma = 0.5E6
            amplitude = self.T0
            y0 = int(self.ny/2)
            x0 = int(self.nx/2)

            for j in range(self.ny):
                for i in range(self.nx):
                    ypart = (self.dy*(j-y0))
                    xpart = (self.dx*(i-x0))
                    self.Pert[j,i] = amplitude * np.exp( - ( xpart**2 + ypart**2 ) / (2*sigma**2))

            """
            two side by side bombs in x-direction
            """
            # sigma = 0.4E6
            # amplitude = 1.7*self.T0
            # y0 = 0
            # #y0 = int(self.ny/2)
            # x0 = int(4*self.nx/5)

            # for j in range(self.ny):
            #     for i in range(self.nx):
            #         ypart = (self.dy*(j-y0))
            #         xpart = (self.dx*(i-x0))
            #         self.Pert[j,i] = amplitude * np.exp( - ( xpart**2 + ypart**2 ) / (2*sigma**2))

            # sigma2 = 0.4E6
            # amplitude2 = 1.7*self.T0
            # y0_2 = 0
            # #y0_2 = int(self.ny/2)
            # x0_2 = int(self.nx/5)

            # for j in range(self.ny):
            #     for i in range(self.nx):
            #         ypart_2 = (self.dy*(j-y0_2))
            #         xpart_2 = (self.dx*(i-x0_2))
            #         self.Pert2[j,i] = amplitude2 * np.exp( - ( xpart_2**2 + ypart_2**2 ) / (2*sigma2**2))



        #Here we add the pertubation to the temperature matrix
        self.T = self.T + self.Pert + self.Pert2

        #Second option, adding the pressure late
        #self.P = self.P0*(self.T/self.T0)**(1/nabla)


       #The other quantites follows from the T-matrix and the P-matrix
        self.e = self.P/(gamma - 1)
        self.rho = self.e * my*mu*(gamma - 1)/(kB*self.T)

        #from here, we have the first enviroment with all the values of the
        #primary variables, aswell as the secondary, including all the boundaries.

    def timestep_fast(self):

        """
        timestep function
        """

        p = 0.1

        rel_rho = np.nanmax(np.abs(self.drhodt / self.rho))
        rel_e = np.nanmax(np.abs(self.dedt / self.e))

        #takes into account that velocity can be zero

        y_u, x_u  = np.where(abs(self.rho_u) > 1e-5)
        y_w, x_w   = np.where(abs(self.rho_w) > 1e-5)

        if len(y_u) > 0:
            rel_u = np.nanmax(abs(self.drho_udt[y_u, x_u] / self.rho_u[y_u, x_u]))
        else:
            rel_u = 0

        if len(y_w) > 0:
            rel_w = np.nanmax(abs(self.drho_wdt[y_w,x_w]/self.rho_w[y_w,x_w]))
        else:
            rel_w = 0

        #Making sure we do not skip a grid point
        rel_x = np.nanmax(np.abs(self.rho_u / (self.rho * self.dx)))
        rel_y = np.nanmax(np.abs(self.rho_w / (self.rho * self.dy)))

        #
        rel_x = np.nanmax(np.abs(self.u / (self.rho * self.dx)))
        rel_y = np.nanmax(np.abs(self.rho_w / (self.rho * self.dy)))

        #Here we determine the largest relative change per time step for any
        #of the quantities at any point on the grid
        max_num = np.array([rel_rho, rel_e, rel_u, rel_w, rel_x, rel_y]).max()

        if max_num == 0:
            dt = p
        else:
            dt = p / max_num

        if dt > p:
            dt = p

        elif dt < 1e-2:     #to ensure dt does not get too small
            dt = 1e-2

        else:
            dt = dt

        return dt


    def boundary_conditions(self):



        """
        boundary conditions for energy, density and velocity

        Here we take care of the vertical boundaries for the n+1 timestep

        """

        my, mu, kB, g, gamma = self.my, self.mu, self.kB, self.g, self.gamma



        alpha_0 = mu*my*g / (kB * self.T[0,:])         #helping variable

        alpha_minus1 = mu*my*g / (kB * self.T[-1,:])   #helping variable


        #internal energy per volume
        self.e[0,:] = (4*self.e[1,:] - self.e[2,:])/(3 - 2*self.dy*alpha_0)
        self.e[-1,:] = (4*self.e[-2,:] - self.e[-3,:])/(3 + 2*self.dy*alpha_minus1)

        ##density

        self.rho[0,:] = self.e[0,:] * my*mu*(gamma-1)/(kB*self.T[0,:])
        self.rho[-1,:] = self.e[-1,:] * my*mu*(gamma-1)/(kB*self.T[-1,:])

        ##velocity

        #u
        self.u[0,:] = (4/3)* self.u[1,:] - (1/3)* self.u[2,:]
        self.u[-1,:] = (4/3) * self.u[-2,:] - (1/3) * self.u[-3,:]

        #w
        self.w[0,:] = 0
        self.w[-1,:] = 0


    def central_x(self,func):


        """
        central difference scheme in x-direction
        """

        right_matrix = np.roll(func,-1, axis= 1) #roll function to i+1, including boundaries
        left_matrix  = np.roll(func,1, axis= 1) #roll function to i-1, including boundaries

        #central_x scheme, now excluding boundary conditions
        return_matrix = (right_matrix[1:-2,:] - left_matrix[1:-2,:])/(2*self.dx)
        return return_matrix


    def central_y(self,func):


        """
        central difference scheme in y-direction
        """

        right_matrix = np.roll(func,-1, axis= 0) #roll function to j+1, including boundaries
        left_matrix  = np.roll(func,1, axis= 0) #roll function to j-1, including boundaries

        #central_y scheme, now excluding boundary conditions

        return_matrix = (right_matrix[1:-2,:] - left_matrix[1:-2,:])/(2*self.dy)
        return return_matrix



    def upwind_x(self,func,u):

        """
        upwind difference scheme in x-direction
        """

        #when velocity is postiive
        left_matrix  = np.roll(func,1, axis= 1) #roll function to i-1, including boundaries

        return_matrix = (func[1:-2,:] - left_matrix[1:-2,:])/self.dx

        #Did we have any negative velocity components?
        y, x = np.where(u[1:-2,:] < 0)

        if len(y) != 0:
            right_matrix = np.roll(func,-1, axis= 1) #roll function to i+1, including boundaries

            return_matrix[y,x] =  (right_matrix[1:-2,:][y,x] - func[1:-2,:][y,x])/self.dx #fixing the negative values using the right schematic


        return return_matrix


    def upwind_y(self,func,u):


        """
        upwind difference scheme in y-direction
        """
        left_matrix = np.roll(func,1, axis = 0) #roll function to j-1, including boundaries

        return_matrix = (func[1:-2,:] - left_matrix[1:-2,:])/self.dy



        ##Did we have any negative velocity components?

        y, x = np.where(u[1:-2,:] < 0)

        if len(y) != 0:

            right_matrix = np.roll(func,-1, axis= 0) #roll function to j+1, including boundaries

            return_matrix[y,x] = (right_matrix[1:-2,:][y,x] - func[1:-2,:][y,x])/self.dy #fixing the negative values using the right schematic


        return return_matrix



    def hydro_solver(self):

        """

        hydrodynamic equations solver

        """

        rho, u, w, e = self.rho, self.u, self.w, self.e

        self.rho_u = rho * u
        self.rho_w = rho * w

        #Continuity equation
        self.drhodt[1:-2,:] = -rho[1:-2,:] *(self.central_x(u) + self.central_y(w)) - u[1:-2,:]*self.upwind_x(rho,u)-w[1:-2,:] * self.upwind_y(rho,w)


        #Momentum equation
        self.drho_udt[1:-2,:] = - self.rho_u[1:-2,:] * (self.upwind_x(u,u) + self.upwind_y(w,u)) - u[1:-2,:]*self.upwind_x(self.rho_u,u) - w[1:-2,:]*self.upwind_y(self.rho_u,w) - self.central_x(self.P)
            #As we have defined y = 0 to be at the bottom of the box, positive y direction is outwards towards the surface. g is pointing in the opposite direction, hence we must have g negative.
        self.drho_wdt[1:-2,:] = - w[1:-2,:] * (self.upwind_y(self.rho_w,w) + self.upwind_x(self.rho_u,w)) - self.rho_w[1:-2,:] * self.upwind_y(w,w)  - self.rho_u[1:-2,:] * self.upwind_x(w,u) - self.central_y(self.P) - rho[1:-2,:]*self.g


        #Energy equation
        self.dedt[1:-2,:] = -e[1:-2,:]* (self.central_x(u) + self.central_y(w)) - u[1:-2,:] * self.upwind_x(e,u) - w[1:-2,:]*self.upwind_y(e,w) - (self.P[1:-2,:] * (self.central_x(u) + self.central_y(w)))
  

        #calculate timestep
        dt = self.timestep_fast()

        #Updating variables for the next time frame
        self.rho[:] = self.rho + self.drhodt*dt
        self.e[:]   = self.e + self.dedt * dt
        self.u[:]   = (self.rho_u + self.drho_udt *dt)/self.rho
        self.w[:]   = (self.rho_w + self.drho_wdt *dt)/self.rho

        #initiate boundary conditions for primary variables
        self.boundary_conditions()

        #Updating secondary variables
        self.T[:] = self.e * (self.gamma-1)*self.mu*self.my/(self.kB*self.rho)
        self.P[:] = (self.gamma-1)*self.e

        return dt
