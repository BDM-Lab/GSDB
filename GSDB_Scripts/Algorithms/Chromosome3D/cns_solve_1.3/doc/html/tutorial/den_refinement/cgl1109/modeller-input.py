from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='hp33_promals3d.ali',
              knowns='1VGYA', sequence='hp33',
              assess_methods=(assess.DOPE, assess.GA341))
a.very_fast()
a.starting_model = 1
a.ending_model = 1
a.make()