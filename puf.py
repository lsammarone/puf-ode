import math
from scipy.integrate import ode
import matplotlib.pyplot as plt
import os

class Pulse:

    def __init__(self,width, height, delay=0.0):
        self.delay = delay
        self.pulse_width = width
        self.height = height

    def compute(self,time):
        delta = time - self.delay
        if delta >= 0 and delta < self.pulse_width:
            return self.height
        else:
            return 0.0

    def copy(self):
        return Pulse(self.pulse_width, self.height, self.delay)

class Zero:

	def __init__(self):
		pass

	def compute(self,time):
		return 0.0;

class ODEPuf3:

    def __init__(self, ks):
        assert(len(self.var_names) == self.number_diffeqs)

        self.ks = {}
        for vn in self.var_names:
            assert(vn in ks)
            self.ks[vn] = ks[vn]


        self._state = {}
        for vn in self.var_names:
            self._state[vn] = 0.0

        self._forcing_functions = {};
        for vn in self.var_names:
            self._forcing_functions[vn] = Zero()


    def state(self):
        return dict(self._state)

    def set_state(self,state):
        for variable in self.var_names:
            self._state[variable] = state[variable]

    @property
    def number_diffeqs(self):
        return self.number_cells*2;

    @property
    def number_cells(self):
        return 3;

    @property
    def var_names(self):
        return ["x0","v0","x1","v1","x2","v2"]

    def set_to_zero(self):
        for variable in self.var_names:
            self._state[variable] = 0.0

    def _derivative(self,t,state):
		#FIXME, return array of derivatives given the current state of the system.
        # x_0
        ddts = dict(map(lambda variable: (variable,0.0), self.var_names))
        ddts["x0"] = self.ks["x0"]*state["x0"] + self._forcing_functions["x0"].compute(t)
        #FIXME: return actual correct derivatives
        return ddts

    def _to_array(self,value_dict):
        return list(map(lambda name: value_dict[name],self.var_names))

    def _from_array(self,value_list):
        return dict(zip(self.var_names, value_list))

    # simulates the system, using the state variables in `self._state` as the initial condition
    # The simulation runs for `sim_time` units of time and returns a dictionary of computed time
    # and value samples
    def simulate(self, sim_time, samples_per_step=100):
        def dt_func(t,vs):
            result = self._derivative(t, self._from_array(vs))
            return self._to_array(result)

        # compute the time between samples and the total number of samples
        dt = 1.0/samples_per_step
        n = math.ceil(sim_time/dt)

        # setup differential equation solver
        solver = ode(dt_func).set_integrator('zvode', method='bdf')

        # set initial condition
        x0 = self._to_array(self.state())
        solver.set_initial_value(x0)

        # record the state at each sample step. The solver may compute the step
        # at intermediate times to get a good simulation.
        info = {"time":[], "values":[]}
        while solver.successful() and solver.t < sim_time:
            info["time"].append(solver.t)
            info["values"].append(solver.y)
            solver.integrate(solver.t + dt)


        return info


    # saves the variable trajectories from the simulation to a series of files in plots/
    # The `prefix` argument sets the prefix for each plot. The suffix of the plot
    # is the state variable the plot shows.
    def save_plots(self,prefix,info):
        directory = "plots"
        if not os.path.exists(directory):
            os.mkdir(directory)

        # translate the data in `info` to a format that can be easily plotted
        times = info["time"]
        variables = {}
        for idx,variable in enumerate(self.var_names):
            variables[variable] = []
            for vect in info["values"]:
                variables[variable].append(vect[idx])

        # plot each variable's trajectory and save it to a file
        for name,values in variables.items():
            plt.plot(times,values)
            plt.savefig("%s/%s_%s.png" % (directory,prefix,name))
            plt.clf()


    def write_challenge(self,pulses):
        # make sure the inputs make sense
        if not (len(pulses.keys()) == self.number_diffeqs):
            raise Exception("expected <%d> pulses, since there are <%d> diffeqs" \
                            % (self.number_diffeqs))
        if not (all(map(lambda p: isinstance(p, Pulse), pulses.values()))):
            raise Exception("the challenge must be a list of Pulse objects")

        # update forcing function array and set delays so that we align signals
        for variable in self.var_names:
            if variable in pulses:
                self._forcing_functions[variable] = pulses[variable].copy()
                self._forcing_functions[variable].delay = 0
            else:
                self._forcing_functions[variable] = Zero()

        # set the state variable values to zero
        self.set_to_zero();

        # simulate for 10 simulation time units
        #FIXME.. how do i simulate
        info = self.simulate(10);

        # set the current state of the state variabes in the system
        # to the last value of each state variable in the simulation
        self.set_state(self._from_array(info["values"][-1]))
        return self.state,info

    def get_response(self):
        #FIXME.. how do i get a response
        # simulate for ten simulation time units and return nothing
        info = self.simulate(10);
        return None,info


def test():
    print("making 3-cell puf....")
    # instantiate the puf with a K parameter for each state variable.
    # the constructor will check to make sure every variable has a k parameter.
    puf3 = ODEPuf3({"x0":0.1,"x1": 0.2, "x2": 0.3, "v0": 1.0, "v1": 0.24, "v2":0.01})

    # instantiate state variables to zero
    puf3.set_to_zero()

    # write a challenge to the puf
    # each puf is labelled with the variable it is provided to
    state,info = puf3.write_challenge({ \
                                        "x0":Pulse(width=1.0,height=1.0), \
                                        "v0":Pulse(width=1.0,height=2.0), \
                                        "x1":Pulse(width=1.0,height=3.0), \
                                        "v1":Pulse(width=2.0,height=1.0), \
                                        "x2":Pulse(width=3.0,height=1.0), \
                                        "v2":Pulse(width=0.5,height=1.0)})

    # save the variable trajectories for the challenge phase
    # all plots are written to the "plots/" directory
    puf3.save_plots("challenge",info);

    # collect the response to the challenge
    response,info = puf3.get_response()

    # save the variable trajectories for the response phase
    # all plots are written to the "plots/" directory
    puf3.save_plots("response",info);

test()
