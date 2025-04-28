import pygpc
from pygpc.AbstractModel import AbstractModel
import numpy as np
from collections import OrderedDict

# Following this tutorial:
# https://github.com/pygpc-polynomial-chaos/pygpc/blob/master/tutorials/Tutorial_pygpc.ipynb

class Ishigami(AbstractModel):
    def validate(self):
        pass

    def simulate(self, process_id=None, matlab_engine=None):
        if self.p["x1"] is not np.ndarray:
            self.p["x1"] = np.array(self.p["x1"])

        if self.p["x2"] is not np.ndarray:
            self.p["x2"] = np.array(self.p["x2"])

        if self.p["x3"] is not np.ndarray:
            self.p["x3"] = np.array(self.p["x3"])

        if self.p["a"] is not np.ndarray:
            self.p["a"] = np.array(self.p["a"])

        if self.p["b"] is not np.ndarray:
            self.p["b"] = np.array(self.p["b"])

        y = (np.sin(self.p["x1"].flatten()) + 
             self.p["a"].flatten() * np.sin(self.p["x2"].flatten()) ** 2
             + self.p["b"].flatten() * self.p["x3"].flatten() ** 4 *
             np.sin(self.p["x1"].flatten()))

        if type(y) is not np.ndarray:
            y = np.array([y])

        y_out = y[:, np.newaxis]

        return y_out

# Secondly, a script is needed to define and run the gPC. In this we define:
#  - the system under investigation (model)
#  - random parameters (problem)
#  - options of the gPC algorithm (algorithm)
#
# The struture of pygpc is:
#    algorithm(problem(model))
#
# This allows you to:
#  - easily replace the model while keeping the problem and the algorithm
#  - define multiple problems using the same model and algorithm (e.g. when different sets of parameters are assumed to be uncertain)  <-- This might be of interest: (a) Uncertainty due to uncertainty of conductivities, (b) Uncertainty due to uncertainty of coil positions
#  - easily replace the algorithm while keeping the problem and the model

model = Ishigami()

parameters = OrderedDict()
parameters["x1"] = pygpc.Beta(pdf_shape = [1, 1], pdf_limits=[-np.pi, np.pi])
parameters["x2"] = pygpc.Beta(pdf_shape = [1, 1], pdf_limits=[-np.pi, np.pi])
parameters["x3"] = 0.
parameters["a"]  = 7.
parameters["b"]  = 0.1

problem = pygpc.Problem(model, parameters)


options = {}
options["order_start"] = 5
options["order_end"] = 20
options["solver"] = "LarsLasso"
options["interaction_order"] = 2
options["order_max_norm"] = 0.7
options["n_cpu"] = 0
options["adaptive_sampling"] = True
options["gradient_enhanced"] = True
options["fn_results"] = "tmp/mygpc" 
options["error_type"] = "nrmsd"
options["eps"] = 0.001

algorithm = pygpc.RegAdaptive(problem=problem, options=options)

session = pygpc.Session(algorithm)

session, coeffs, results = session.run()

mean = session.gpc[0].get_mean(coeffs)
print(mean)
std  = session.gpc[0].get_std(coeffs)
print(std)
