import pygpc
from pygpc.AbstractModel import AbstractModel
import numpy as np
from collections import OrderedDict


class DistanceModel(AbstractModel):
    def __init__(self):
        super().__init__()
        self.how_often_evaluated = 0

    def validate(self):
        pass

    def simulate(self, process_id=None, matlab_engine=None):
        if self.p["x1"] is not np.ndarray:
            self.p["x1"] = np.array(self.p["x1"]).flatten()

        if self.p["x2"] is not np.ndarray:
            self.p["x2"] = np.array(self.p["x2"]).flatten()

        if self.p["y1"] is not np.ndarray:
            self.p["y1"] = np.array(self.p["y1"]).flatten()

        if self.p["y2"] is not np.ndarray:
            self.p["y2"] = np.array(self.p["y2"]).flatten()

        # The entries of self.p are arrays with shape (n,)

        print(f"\n\n\n\n\n\np: {self.p}")
        print(f"shape of x1: {self.p['x1'].shape}")
        print(f"shape of x2: {self.p['x2'].shape}")
        print(f"shape of y1: {self.p['y1'].shape}")
        print(f"shape of y2: {self.p['y2'].shape}")

        n_eval = len(self.p["x1"].flatten())
        d = np.zeros((n_eval,)) * np.nan
        for i in range(n_eval):
            d1 = (self.p["x1"].flatten()[i] - self.p["y1"].flatten()[i])**2
            d2 = (self.p["x2"].flatten()[i] - self.p["y2"].flatten()[i])**2
            d[i] = np.sqrt(d1 + d2)
            self.how_often_evaluated += 1

        #print(f"\nshape of d: {d.shape}")

        d_out = d[:, np.newaxis]
        print(f"\nshape of d_out: {d_out.shape}")

        return d_out
    


model = DistanceModel()

parameters = OrderedDict()
parameters["x1"] = pygpc.Norm(pdf_shape=[0, 0.5])
parameters["x2"] = pygpc.Norm(pdf_shape=[0, 0.5]) # I.e. \vec{x} is somewhere around the origin (0,0)
parameters["y1"] = pygpc.Norm(pdf_shape=[0, 0.5])
parameters["y2"] = pygpc.Norm(pdf_shape=[1, 0.5]) # I.e. \vec{y} is somewhere around (0, 1)  --- approx. 1 unit away

problem = pygpc.Problem(model, parameters)


options = {}
options["order_start"] = 5
options["order_end"] = 20
options["solver"] = "LarsLasso"
options["interaction_order"] = 2
options["order_max_norm"] = 0.7
options["n_cpu"] = 0
options["adaptive_sampling"] = True
options["gradient_enhanced"] = False
options["fn_results"] = "tmp/mygpc" 
options["error_type"] = "loocv"
options["error_norm"] = "absolute"
options["eps"] = 0.1
options["GPU"] = True

algorithm = pygpc.RegAdaptive(problem=problem, options=options)

session = pygpc.Session(algorithm)

import warnings
warnings.filterwarnings("ignore")

session, coeffs, results = session.run()

print(f"How often was the model evaluated: {model.how_often_evaluated} times")

mean = session.gpc[0].get_mean(coeffs)
print(f"mean = {mean}")
std  = session.gpc[0].get_std(coeffs)
print(f"std = {std}")
