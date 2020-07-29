import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches
import random


def fluxGraph(nt, dt, peopleConfig, title, haveModule, haveHVAC):
	# things set constant for this experiment
	nx = 120
	ny = 40
	nit = 50
	c = 1
	dx = 2 / (nx - 1)
	dy = 2 / (ny - 1)
	x = numpy.linspace(0, int(nx/20), nx)
	y = numpy.linspace(0, int(ny/20), ny)
	X, Y = numpy.meshgrid(x, y)
	rho = 3
	nu = .1
	fanSpeed = 1.59
	HVACSpeed = 3
	addS = "withModule"
	if not haveModule:
		addS = "withoutModule"
		fanSpeed = 0
	addX = "withHVAC"
	if not haveHVAC:
		addX = "withoutHVAC"
		HVACSpeed = 0



	u = numpy.zeros((ny, nx))
	v = numpy.zeros((ny, nx))
	p = numpy.zeros((ny, nx))
	b = numpy.zeros((ny, nx))


	u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, peopleConfig, fanSpeed, HVACSpeed, nx,ny)
	fig = pyplot.figure(figsize=(13,6), dpi=100)

	u_sums = [sum(u[:,i]) for i in range(nx)]

	pyplot.plot(x, u_sums)

	pyplot.savefig('_fluxStuff_' + addS + '_' + addX + '_' + str(nt*dt) + '_' + title + '.png' )
	pyploy.close()

# plotType = ‘stream’ or ‘vector’
# nt is the number of timesteps
# dt is the timestep length
# peopleConfig is array, storing 1 or 0 based on where people are
#	eg [0,1,1,0] means two people sitting
# title is the title that the thing is saved as
def ventExperiment(plotType, nt, dt, peopleConfig, title, haveModule, haveHVAC):
	# things set constant for this experiment
	nx = 120
	ny = 40
	nit = 50
	c = 1
	dx = 2 / (nx - 1)
	dy = 2 / (ny - 1)
	x = numpy.linspace(0, int(nx/20), nx)
	y = numpy.linspace(0, int(ny/20), ny)
	X, Y = numpy.meshgrid(x, y)
	rho = 3
	nu = .1
	fanSpeed = 1.59
	HVACSpeed = 3
	addS = "withModule"
	if not haveModule:
		addS = "withoutModule"
		fanSpeed = 0
	addX = "withHVAC"
	if not haveHVAC:
		addX = "withoutHVAC"
		HVACSpeed = 0

	u = numpy.zeros((ny, nx))
	v = numpy.zeros((ny, nx))
	p = numpy.zeros((ny, nx))
	b = numpy.zeros((ny, nx))


	u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, peopleConfig, fanSpeed, HVACSpeed, nx,ny)
	fig = pyplot.figure(figsize=(13,6), dpi=100)

	table1 = patches.Rectangle((1.0,0.4), 1, 0.1)
	table2 = patches.Rectangle((4.0,0.4), 1, 0.1)
	if haveModule:
		module1 = patches.Rectangle((1.5,0.5), 0.1, 0.1)
		module2 = patches.Rectangle((4.5,0.5), 0.1, 0.1)
		people = [table1,table2, module1, module2]
	else:
		people = [table1,table2]
	for i in range(len(peopleConfig)):
		if peopleConfig[i] == 1:
			if i == 0:
				people.append(patches.Rectangle((0.6,0.4), 0.2, 0.3))
			if i == 1:
				people.append(patches.Rectangle((2.2,0.4), 0.2, 0.3))
			if i == 2:
				people.append(patches.Rectangle((3.6,0.4), 0.2, 0.3))
			if i == 3:
				people.append(patches.Rectangle((5.2,0.4), 0.2, 0.3))

	print("lem", str(len(people)))

	if plotType == 'stream':
		fig = pyplot.figure(figsize=(13,6), dpi=100)
		pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)
		pyplot.colorbar()
		pyplot.contour(X, Y, p, cmap=cm.viridis)
		pyplot.streamplot(X, Y, u, v)
		pyplot.xlabel('X')
		pyplot.ylabel('Y');
	elif plotType == 'vector':
		fig = pyplot.figure(figsize=(13,6), dpi=100)
		# plotting the pressure field as a contour
		pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)
		pyplot.colorbar()
		# plotting the pressure field outlines
		pyplot.contour(X, Y, p, cmap=cm.viridis)
		# plotting velocity field
		pyplot.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2])
		pyplot.xlabel('X')
		pyplot.ylabel('Y')

	ax = pyplot.gca()
	for p in range(len(people)):
		ax.add_patch(people[p])
		pyplot.title("sds")
		pyplot.title('t = ' + str(nt*dt) + 's')

	pyplot.savefig(plotType + '_vent_' + addS + '_' + addX + '_' + str(nt*dt) + '_' + title + '.png' )
	pyplot.close()






def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, peopleConfig, fanSpeed, HVACSpeed,nx, ny):
    un = numpy.empty_like(u)
    vn = numpy.empty_like(v)
    b = numpy.zeros((ny, nx))

    for n in range(nt):
        un = u.copy()
        vn = v.copy()

        b = build_up_b(b, rho, dt, u, v, dx, dy)
        p = pressure_poisson(p, dx, dy, b)

        u[1:-1, 1:-1] = (un[1:-1, 1:-1]-
                         un[1:-1, 1:-1] * dt / dx *
                        (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * dt / dy *
                        (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                         nu * (dt / dx**2 *
                        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                         dt / dy**2 *
                        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

        v[1:-1,1:-1] = (vn[1:-1, 1:-1] -
                        un[1:-1, 1:-1] * dt / dx *
                       (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                        vn[1:-1, 1:-1] * dt / dy *
                       (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                        nu * (dt / dx**2 *
                       (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                        dt / dy**2 * (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))
        AC_bot = 2
        AC_top = 4
        u[0:AC_bot, 0]  = 0 #influx of air
        u[AC_bot:AC_top, 0]  = HVACSpeed #influx of air
        u[AC_top:-1, 0]  = 0 #influx of air

        #AC_bot = 2  + 30
        #AC_top = 5 + 30

        u[0:AC_bot, -1]  = 0 #influx of air
        u[AC_bot:AC_top, -1]  = HVACSpeed #influx of air
        u[AC_top:-1, -1]  = 0 #influx of air

        #table 1
        v[8:10, 20:40] = 0
        u[8:10, 20:40] = 0

        #table 2
        v[8:10, 80:100] = 0
        u[8:10, 80:100] = 0

        #rect3 = patches.Rectangle((2.2,0.4), 0.2, 0.3)
        #rect4 = patches.Rectangle((3.6,0.4), 0.2, 0.3)
        v[11,29:31]=fanSpeed
        v[11,89:91] =fanSpeed


        for i in range(len(peopleConfig)):
            if peopleConfig[i] == 1:
                if i == 0:
                    v[8:14,12:16] = 0
                    u[8:14,12:16] = 0
                    u[12, 18] = 2.2*numpy.sin(n * 2*3.1415/5000 + random.random())
                if i == 1:
                    v[8:14,44:48] = 0
                    u[8:14,44:48] = 0
                    #breathe
                    u[12, 42] = 2.2*numpy.sin(n * 2*3.1415/5000 + random.random() )
                if i == 2:
                    v[8:14,72:76] = 0
                    u[8:14,72:76] = 0
                    #breathe
                    u[12, 78] = 2.2*numpy.sin(n * 2*3.1415/5000 + random.random())
                if i == 3:
                    v[8:14,104:108] = 0
                    u[8:14,104:108] = 0
                    #breathe
                    u[12, 102] = 2.2*numpy.sin(n * 2*3.1415/5000 + random.random() )

        #The REST OF THE BOUNDARY CONDITIONS
        u[-1, :] = 0
        u[0, :]  = 0
        v[0, :]  = 0
        v[-1, :] = 0
        v[:, 0]  = 0
        v[:, -1] = 0

    return u, v, p


def pressure_poisson(p, dx, dy, b):
    pn = numpy.empty_like(p)
    pn = p.copy()

    nit = 50
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
                          (2 * (dx**2 + dy**2)) -
                          dx**2 * dy**2 / (2 * (dx**2 + dy**2)) *
                          b[1:-1,1:-1])

        p[:, -1] = p[:, -2] # dp/dx = 0 at x = 2
        p[0, :] = p[1, :]   # dp/dy = 0 at y = 0
        p[:, 0] = p[:, 1]   # dp/dx = 0 at x = 0
        p[-1, :] = 0        # p = 0 at y = 2

    return p

def build_up_b(b, rho, dt, u, v, dx, dy):

    b[1:-1, 1:-1] = (rho * (1 / dt *
                    ((u[1:-1, 2:] - u[1:-1, 0:-2]) /
                     (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                    ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 -
                      2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                           (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))-
                          ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))

    return b


ventExperiment('stream', 1000, 0.001, [0,1,1,0], 'test1.1', True, True)
ventExperiment('stream', 1000, 0.001, [0,1,0,1], 'test1.2', True, True)
ventExperiment('stream', 1000, 0.001, [1,0,0,1], 'test1.3', True, True)
ventExperiment('stream', 1000, 0.001, [1,0,1,0], 'test1.4', True, True)
ventExperiment('stream', 1000, 0.001, [0,1,1,0], 'test1.5', False, False)
ventExperiment('stream', 1000, 0.001, [1,0,0,1], 'test1.6', False, False)
ventExperiment('stream', 1000, 0.001, [0,1,1,0], 'test1.7', True, False)
fluxGraph(1000, 0.001, [0,1,1,0], 'test1.7', True, False)
fluxGraph(1000, 0.001, [0,1,1,0], 'test1.1', True, True)
