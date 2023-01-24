"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import os

src = os.path.dirname(os.path.abspath(__file__))
codebase = os.path.abspath(os.path.join(src, os.pardir))
projectroot = os.path.abspath(os.path.join(codebase, os.pardir))
figures = os.path.abspath(os.path.join(projectroot, 'figures'))
figures_export = os.path.abspath(os.path.join(codebase, 'export figs'))

data = os.path.abspath(os.path.join(codebase, 'data'))
vault = os.path.abspath(os.path.join(data, 'vault'))
tests = os.path.abspath(os.path.join(codebase, 'tests'))
fixtures = os.path.abspath(os.path.join(tests, 'fixtures'))

eqexp = os.path.abspath(os.path.join(data, 'equal expectancy analysis'))
predictions = os.path.abspath(os.path.join(eqexp, 'predicted scores'))
model_results = os.path.abspath(os.path.join(
    eqexp, 'model results from predictied'))
