"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2023, B.Sz.
"""

import os

src = os.path.dirname(os.path.abspath(__file__))
codebase = os.path.abspath(os.path.join(src, os.pardir))
projectroot = os.path.abspath(os.path.join(codebase, os.pardir))
figures = os.path.abspath(os.path.join(projectroot, 'figures'))
#figures_export = os.path.abspath(os.path.join(codebase, 'export figs'))
figures_export = os.path.abspath(os.path.join(codebase, 'export figs', 'private', 'tmp'))
data = os.path.abspath(os.path.join(codebase, 'data'))
