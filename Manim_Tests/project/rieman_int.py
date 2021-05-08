
from manim import *

from stareg.bspline import Bspline
import numpy as np


class RiemannRectanglesAnimation(GraphScene):
    def __init__(self, **kwargs):
            GraphScene.__init__(
                self,
                y_min=0,
                y_max=1,
                x_min=-0.1,
                x_max=1.1,
                **kwargs
            )
            
    arguments = {
        "y_max": 1,
        "y_axis_height": 1,
        "x_min":0,
        "x_max":1,
        "init_dx":0.05,
    }
    def construct(self):
        self.setup_axes()
        def func(x):
            return 0.1 * (x + 3-5) * (x - 3-5) * (x-5) + 5

        def func2(x):
            BS = Bspline()
            B, k = BS.basismatrix(np.linspace(0,1,100)).values()
            return float(BS.basisfunction(np.asarray(x), k, 5, 3))

        graph=self.get_graph(func2,x_min=0.0,x_max=1.)
        kwargs = {
            "x_min" : 0,
            "x_max" : 1,
            "fill_opacity" : 0.75,
            "stroke_width" : 0.25,
        }
        flat_rectangles = self.get_riemann_rectangles(
                                self.get_graph(lambda x : 0),
                                dx=0.05,
                                start_color=invert_color(PURPLE),
                                end_color=invert_color(ORANGE),
                                **kwargs
        )
        riemann_rectangles_list = self.get_riemann_rectangles_list(
                                graph,
                                2,
                                max_dx=0.05,
                                power_base=2,
                                start_color=PURPLE,
                                end_color=ORANGE,
                                 **kwargs
        )
        self.add(graph)
        # Show Riemann rectangles
        self.play(ReplacementTransform(flat_rectangles,riemann_rectangles_list[0]))
        self.wait()
        for r in range(1,len(riemann_rectangles_list)):
            self.transform_between_riemann_rects(
                    riemann_rectangles_list[r-1],
                    riemann_rectangles_list[r],
                    replace_mobject_with_target_in_scene = True,
                )
        self.wait()