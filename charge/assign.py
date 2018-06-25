import os
import random

import networkx as nx

from charge.babel import convert_from, IOType
from charge.collectors import AtomicMeanCollector, AtomicHistogramCollector
from charge.repository import Repository
from charge.solvers import SimpleSolver, ILPSolver, DPSolver, CDPSolver


class AtomicSimpleCharger(AtomicMeanCollector, SimpleSolver):
    pass


class AtomicILPCharger(AtomicHistogramCollector, ILPSolver):
    pass


class AtomicDPCharger(AtomicHistogramCollector, DPSolver):
    pass


class AtomicCDPCharger(AtomicHistogramCollector, CDPSolver):
    pass
