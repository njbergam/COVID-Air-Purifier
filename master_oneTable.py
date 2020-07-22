import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches



# plotType = ‘stream’ or ‘vector’
# nt is the number of timesteps
# dt is the timestep length
# peopleConfig is array, storing 1 or 0 based on where people are
#	eg [0,1] means two people sitting
# title is the title that the thing is saved as
def oneTable(plotType, nt, dt, peopleConfig, title):
	# things set constant for this experiment
	nx = 100
	ny = 100
	nit = 50
	c = 1
	dx = 2 / (nx - 1)
	dy = 2 / (ny - 1)
	x = numpy.linspace(0, int(nx/25), nx)
	y = numpy.linspace(0, int(ny/25), ny)
	X, Y = numpy.meshgrid(x, y)
	rho = 3
	nu = .1
	fanSpeed = 1.59
	breath  = 0

	u = numpy.zeros((ny, nx))
	v = numpy.zeros((ny, nx))
	p = numpy.zeros((ny, nx))
	b = numpy.zeros((ny, nx))


	u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, peopleConfig, fanSpeed, nx,ny)
	fig = pyplot.figure(figsize=(9,8), dpi=100)

	table = patches.Rectangle((1.0,0.5), 2, 0.1)
	module = patches.Rectangle((1.9,0.6), 0.2, 0.1)

	people = [table,module]
	for i in range(len(peopleConfig)):
		if peopleConfig[i] == 1:
			if i == 0:
				people.append(patches.Rectangle((0.6,0.7), 0.2, 0.3))
			elif i == 1:
				people.append(patches.Rectangle((3.2,0.7), 0.2, 0.3))


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
	pyplot.savefig(plotType + '_vent_' + str(nt*dt) + '_' + title + '.png' )
	pyplot.show()
	pyplot.close()



def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, peopleConfig, fanSpeed, nx, ny):
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
                        dt / dy**2 *
                       (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))


        #rect3 = patches.Rectangle((2.2,0.4), 0.2, 0.3)
        #rect4 = patches.Rectangle((3.6,0.4), 0.2, 0.3)


        for i in range(len(peopleConfig)):
            if peopleConfig[i] == 1:
                if i == 0:
                    v[8:14,15:20] = 0
                    u[8:14,15:20] = 0
                    u[20, 23] = 2.2
                if i == 1:
                    v[8:14,80:85] = 0
                    u[8:14,80:85] = 0
                    #breathe
                    u[20, 77] = 2.2


        # fans
        v[16,50] = fanSpeed

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


oneTable('stream', 1000, 0.001, [1,1], 'sethDad')
