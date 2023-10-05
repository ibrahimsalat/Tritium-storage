import numpy as np 
import matplotlib.pyplot as plt



#initial_pellet_size = 50 #microns
N = 100
rho_u = 79663.866
rho_uh3 =  45435.68 #0.045/(100**3)
b = 0.67
velocity = 0.5
rho_h = 0.081
p = 100000 #1 bar

#every unit cell has 8 uranium 

def visccosity(T):
    return (8.76*10e-6)*(365.85/(T+72))*((T/293.85)**1.5)

def diffusitivity_helium(T):
    return 1.032*1e-8 * (T**1.74)

def reaction_rate(T):
    return 0.51*((8.314*T)**0.5)*np.exp(-3033/T)#1#10.4*np.exp(-1592/T)


def diffusivity_constant(T, a):
    return a*np.exp(-2300/(8.314*T))#5e-8*np.exp(80542/8.314*T)# #2.11e-8*np.exp(-2300/(8.314*T))##1.2e9*np.exp(25522.4/8.314*T)  #2.11e-6*np.exp(-5820/(T))



def core_radius(initial_pellet_size, shrinking_step):
    core_radius_values =  np.arange(initial_pellet_size,-1*shrinking_step , - shrinking_step )
    return core_radius_values
    
#r_c = core_radius(initial_pellet_size, shrinking_step)
def outer_radius(initial_pellet_size, rho_u, rho_uh3, r_c):
    outer_radius_values = np.zeros(N+1)
    for i in range(N+1):
        outer_radius_values[i] =((((initial_pellet_size*1e-6)**3)-((r_c[i]*1e-6)**3))*
                                 (rho_u/rho_uh3)+((r_c[i]*1e-6)**3))**(1/3)
    return outer_radius_values

def mass_transfer(r_o,d, viscosity, r_c):
    kg = np.zeros(N+1)
    sc = viscosity/d
    for i in range(N+1):
        re = rho_h*velocity*2*r_c[i]*1e-6/viscosity
        kg[i] = (2 + (0.6*(sc**(1/3)) * (re**0.5)))*(d)/(2*r_o[i]*1e-6)
    return kg



def const(kr, kg, r_c, r_o, de):
    a =  np.zeros(N+1)
    b =  np.zeros(N+1)
    e =  np.zeros(N+1)
    const = np.zeros(N+1)
    for i in range(N+1):
        a[i] = (((r_c[i]*1e-6)**2)*kr)/((r_o[i]**2)*kg[i])
        b[i] = ((r_c[i]*1e-6)*kr)/(de)
        e[i] = (((r_c[i]*1e-6)**2)*kr)/(de*r_o[i])
        const[i] = a[i]+b[i]-e[i]
    return const



def core_concentration(c_b, const):
    c_c = np.zeros(N+1)
    for i in range(N+1):
        c_c[i] = ((-const[i] + np.sqrt((const[i]**2)+(4*c_b)))/2)**2
    return c_c

#c_c = core_concentration(c_b)
#print(c_c)

def time_taken(b, c_c, shrinking_step, kr):
    time_diff = np.zeros(N+1)
    for i in range(N+1):
        time_diff[i] = ((shrinking_step* 1e-6)*rho_u)/(b*kr*np.sqrt(c_c[i]))#
    
    time = np.cumsum(time_diff)
    return time



def thiele_modulus(initial_pellet_size , kr, c_c, de):
    m = 2*initial_pellet_size*1e6 *(np.sqrt((1.5*kr*(c_c**(-0.5)))/2*de))
    return m


#----------------------------------------------------------------------------
#Initial Conditions with varying pellet size
T = 500

#calling the functions: 
def simulation(T, initial_pellet_size, p, N, a):
    shrinking_step = initial_pellet_size/N
    c_b =  p/(8.314*T)
    d = diffusitivity_helium(T)
    viscosity = visccosity(T)
    kr = reaction_rate(T)  
    de = diffusivity_constant(T, a)
    r_c = core_radius(initial_pellet_size, shrinking_step)
    r_o = outer_radius(initial_pellet_size, rho_u, rho_uh3, r_c)
    kg = mass_transfer(r_o,d, viscosity, r_c)
    const1 = const(kr,kg, r_c, r_o,de)
    c_c = core_concentration(c_b, const1)
    time = time_taken(0.667, c_c, shrinking_step, kr)
    x = r_c/initial_pellet_size
    relitive_sphere_radius = r_o/(initial_pellet_size*1e-6)
    pressure_change = (1800 - ((3/((4/3)*3.14*((initial_pellet_size*1e-6)**3)*19.1e6))*((4/3)*3.14*((initial_pellet_size*1e-6)**3 - (r_c*1e-6)**3))*rho_u*1.5*8.314*T/0.00048)/100)
    m = thiele_modulus(initial_pellet_size , kr, c_c, de)
    y =  np.log10(c_c)
    return time, x, pressure_change, y, relitive_sphere_radius, m

a = 1.11e-10
initial_pellet_size = 2.5
time, x, pressure_change, y, relitive_sphere_radius, m = simulation(T, initial_pellet_size, p, N,a)
initial_pellet_size2 = 3
time2, x2, pressure_change2, y2, relitive_sphere_radius2, m2 = simulation(T,initial_pellet_size2 , p, N,a)
initial_pellet_size3 = 3.5
time3, x3, pressure_change3, y3, relitive_sphere_radius3, m3 = simulation(T, initial_pellet_size3 , p, N,a)
initial_pellet_size4 = 0.8
time4, x4, pressure_change4, y4, relitive_sphere_radius4, m4 = simulation(T, initial_pellet_size4 , p, N,a)
initial_pellet_size5 = 1
time5, x5, pressure_change5, y5, relitive_sphere_radius5, m5 = simulation(T, initial_pellet_size5 , p, N,a)
initial_pellet_size6 = 1.2
time6, x6, pressure_change6, y6, relitive_sphere_radius6, m6 = simulation(T, initial_pellet_size6 , p, N,a)

# #_____
plt.rcParams.update({'font.size': 18})

plt.plot(x, y, label = str(initial_pellet_size) + " $\mathregular{\mu}m$")
plt.plot(x2, y2, label = str(initial_pellet_size2) + " $\mathregular{\mu}m$")
plt.plot(x3, y3, label = str(initial_pellet_size3) + " $\mathregular{\mu}m$")
#plt.title("Concentration of the Gas on the Core at " + str(T)+ "K")
plt.xlabel(" Relative Core Radius, $\mathregular{r_c/R_0}$", fontsize = "20")
plt.ylabel("Concentration on the Core, Log[$\mathregular{mol/m^{-3}}$]",fontsize = "18")
plt.legend()
plt.show()


plt.plot(time, pressure_change, label = str(initial_pellet_size) + " $\mathregular{\mu}m$")
plt.plot(time2, pressure_change2, label = str(initial_pellet_size2) + " $\mathregular{\mu}m$")
plt.plot(time3, pressure_change3, label = str(initial_pellet_size3) + " $\mathregular{\mu}m$")
plt.ylabel("Pressure (mbar)")
plt.xlabel("Time (s)")
#plt.xlim([0,1500])
plt.legend()
plt.show()




plt.plot(x, m, label = str(initial_pellet_size) + "$\mathregular{\mu}m$")
plt.plot(x2, m2, label = str(initial_pellet_size2) + "$\mathregular{\mu}m$")
plt.plot(x3, m3, label = str(initial_pellet_size3) + "$\mathregular{\mu}m$")
plt.ylabel("Thiele Modulus")
plt.xlabel("Relative Core Radius, $\mathregular{r_c/R_0}$")
#plt.xlim([0,1500])
plt.legend(loc = 'lower center')
plt.show()

plt.xlabel("Time (s)", fontsize = '20')
plt.ylabel("Relative Core Dimension ($\mathregular{r_c/R_0})$", fontsize = '20')
plt.plot(time, x, label = str(initial_pellet_size) + " $\mathregular{\mu}m$")
plt.plot(time2, x2, label = str(initial_pellet_size2) + " $\mathregular{\mu}m$")
plt.plot(time3, x3, label = str(initial_pellet_size3) + " $\mathregular{\mu}m$")
# plt.plot(time4, x4, label = str(initial_pellet_size4) + " $\mathregular{\mu}m$")
# plt.plot(time5, x5, label = str(initial_pellet_size5) + " $\mathregular{\mu}m$")
# plt.plot(time6, x6, label = str(initial_pellet_size6) + " $\mathregular{\mu}m$")

plt.legend()
plt.show()


plt.plot(time, relitive_sphere_radius, label = str(initial_pellet_size) + " $\mathregular{\mu}m$")
plt.plot(time2 , relitive_sphere_radius2, label = str(initial_pellet_size2) + "$\mathregular{\mu}m$")
plt.plot(time3 , relitive_sphere_radius3, label = str(initial_pellet_size3) + " $\mathregular{\mu}m$")
plt.xlabel("Time (s)")
plt.ylabel("Relative Dimension of the Particle")
plt.legend(loc ="lower right")
plt.show()
#----------------------------

#----------------------------------------------------------------------------
initial_pellet_size = 50
a = 2.11e-8
T = 300
time, x, pressure_change, y, relitive_sphere_radius, m = simulation(T, initial_pellet_size, p, N,a)
T2 = 400
time2, x2, pressure_change2, y2, relitive_sphere_radius2, m2 = simulation(T2, initial_pellet_size, p, N,a)
T3 = 500
time3, x3, pressure_change3, y3, relitive_sphere_radius3, m3 = simulation(T3, initial_pellet_size, p, N,a)
T4 = 600
time4, x4, pressure_change4, y4, relitive_sphere_radius4, m4 = simulation(T4, initial_pellet_size, p, N,a)
T5 = 700
time5, x5, pressure_change5, y5, relitive_sphere_radius5, m5 = simulation(T5, initial_pellet_size, p, N,a)


plt.plot(x, y, label = str(T) + "K")
plt.plot(x, y2, label = str(T2) + "K")
plt.plot(x, y3, label = str(T3) + "K")
plt.xlabel(" Relative Core Radius,  $\mathregular{r_c/R_0}$")
plt.ylabel("Concentration on the core, Log[$\mathregular{mol/m^{-3}}$]")
plt.title("Concentration of the Gas on a " +str(initial_pellet_size) + " micron Particle")
plt.legend()
plt.show()


plt.plot(time, pressure_change, label = str(T) + "K")
plt.plot(time2, pressure_change2, label = str(T2) + "K")
plt.plot(time3, pressure_change3, label = str(T3) + "K")
plt.ylabel("moles of Uranium Core")
plt.xlabel("Time (s)")
#plt.xlim([0,1500])
plt.title("Moles of Uranium used through the Reaction")
plt.legend()
plt.show()



plt.plot(x, m, label = str(T) + "K ")
plt.plot(x2, m2, label = str(T2) + "K")
plt.plot(x3, m3, label = str(T3) + "K")
plt.ylabel("Thiele Modulus")
plt.xlabel("Relative Core Radius, $\mathregular{r_c/R_0}$")
#plt.xlim([0,1500])
plt.legend()
plt.show()

plt.xlabel("Time (s)")
plt.ylabel("Relative dimension ($\mathregular{r_c/R_0}$)")
plt.plot(time, x, label = str(T) + "K")
plt.plot(time2, x, label = str(T2) + "K")
plt.plot(time3, x, label = str(T3) + "K")
plt.plot(time4, x, label = str(T4) + "K")
plt.plot(time5, x, label = str(T5) + "K")
#plt.title("Core Radius as a function of Time for " +str(initial_pellet_size) +" micron particles")
plt.legend()
plt.show()
#-------------------------------------------------------
initial_pellet_size = 3
T = 298
a = 2.11e-9
time, x, pressure_change, y, relitive_sphere_radius, m = simulation(T, initial_pellet_size, p, N,a)
a2 = 2.11e-10
time2, x2, pressure_change2, y2, relitive_sphere_radius2, m2 = simulation(T, initial_pellet_size, p, N,a2)
a3 = 2.11e-11
time3, x3, pressure_change3, y3, relitive_sphere_radius3, m3 = simulation(T, initial_pellet_size, p, N,a3)
a4 = (2.11e-11)
time4, x4, pressure_change4, y4, relitive_sphere_radius4, m4 = simulation(T, initial_pellet_size, p, N,a4)
a5 = (2.11e-11)/np.sqrt(2)
time5, x5, pressure_change5, y5, relitive_sphere_radius5, m5 = simulation(T, initial_pellet_size, p, N,a5)
a6 = (2.11e-11)/np.sqrt(3)
time6, x6, pressure_change6, y6, relitive_sphere_radius6, m6 = simulation(T, initial_pellet_size, p, N,a6)

plt.xlabel("Time (s)")
plt.ylabel("Relative dimension ($\mathregular{r_c/R_0}$)")
plt.plot(time, x, label = "A = " +str(a))
plt.plot(time2, x, label = "A = " +str(a2))
plt.plot(time3, x, label = "A = " +str(a3))
# plt.plot(time4, x, label = "A = " +str(a4))
# plt.plot(time5, x, label = "A = " +str(a5))
# plt.plot(time6, x, label = "A = " +str(a6))
# plt.title("Core Radius as a function of Time for " +str(initial_pellet_size) +" micron particles at " +str(T) +" K")
plt.legend()
plt.show()

plt.xlabel("Time (s)")
plt.ylabel("Relative dimension ($\mathregular{r_c/R_0}$)")
plt.plot(time4, x, label = "A = Hydrogen")
plt.plot(time5, x, label = "A = Deuterium" )
plt.plot(time6, x, label = "A = Tritium")
#plt.title("Core Radius as a function of Time for " +str(initial_pellet_size) +" micron particles at " +str(T) +" K")
plt.legend()
plt.show()


#-------------------------------------------------------
    
# plt.scatter(datax[:-1], rate)
# plt.show()
end_time = np.zeros(20)
temp = np.zeros(20)
rate = np.zeros(20)
a = 2.11e-11
for i in range(20):
    time, x, pressure_change, y, relitive_sphere_radius, m = simulation(250+i*25,3, p, N,a)
    end_time[i] = time[-2]
    rate[i] = 1/time[-2]
    temp[i] = 250+i*25
    
end_time2 = np.zeros(20)
temp2 = np.zeros(20)
rate2 = np.zeros(20)
a1 = 2.11e-10
for i in range(20):
    time2, x2, pressure_change2, y2, relitive_sphere_radius2, m2 = simulation(250+i*25,3, p, N,a1)
    end_time2[i] = time2[-2]
    rate2[i] = 1/time2[-2]
    temp2[i] = 250+i*25
    
end_time3 = np.zeros(20)
temp3 = np.zeros(20)
rate3 = np.zeros(20)
a2 = 2.11e-9
for i in range(20):
    time3, x3, pressure_change3, y3, relitive_sphere_radius3, m3 = simulation(250+i*25,3, p, N,a2)
    end_time3[i] = time3[-2]
    rate3[i] = 1/time3[-2]
    temp3[i] = 250+i*25    
    
end_time4 = np.zeros(20)
temp4 = np.zeros(20)
rate4 = np.zeros(20)
a3 = 2.11e-8
for i in range(20):
    time4, x4, pressure_change4, y4, relitive_sphere_radius4, m4 = simulation(250+i*25,3, p, N,a3)
    end_time4[i] = time4[-2]
    rate4[i] = 1/time4[-2]
    temp4[i] = 250+i*25    

#plt.rcParams.update({'font.size': 16})
plt.style.use("bmh")
plt.plot(temp, (rate), label = "A = " + str(a))
# # plt.plot(temp, (rate2), label = "A = "+ str(a1))
# # plt.plot(temp, (rate3), label = "A = "+ str(a2))
# # plt.plot(temp, (rate4), label = "A = "+ str(a3))
plt.xlabel("Temperature (K)", fontsize = "22")
plt.ylabel("Rate ($\mathregular{s^{-1}})$", fontsize = "22")
# #log[($\mathregular{s^{-1}}$)]
# # plt.yscale("log")
plt.legend()
plt.xlim([250,550])
# plt.title("Reaction Rate for 3$\mathregular{\mu}m$ Particles", fontsize = "20")
plt.show()

tem = [25.13108362870971+273,
50.004791013288205+273, 
99.65725167526503+273, 
149.37924127428502+273, 
199.37490107885773+273,
249.6816431657637+273]

et = [ 786.7238374183794,
  675.3890094469282,
  881.1522851008001,
  2107.94900742607,
  2099.9457622975574,
  2023.0872698914459]

plt.scatter(temp, end_time, color='red', linewidths=3, label = "Modelled SC Data")
plt.scatter(tem, et, color='blue', linewidths=3, label = "Experimental Data")
plt.xlim([250,530])
plt.ylim([0,3000])
plt.legend()
plt.xlabel("Temperature (K)")
plt.ylabel("Reaction End Time (s)")
plt.show ()
