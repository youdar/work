#!/usr/bin/env python

from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

from test import Banana


graphviz = GraphvizOutput(output_file='filter_none.png')

with PyCallGraph(output=GraphvizOutput()):
    banana = Banana()
    banana.eat()