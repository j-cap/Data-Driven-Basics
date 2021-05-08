
from manim import *

from stareg.bspline import Bspline
import numpy as np

class BS_Graph(GraphScene):
    def __init__(self, **kwargs):
        GraphScene.__init__(
            self,
            x_min=0,
            x_max=10,
            y_min=0,
            y_max=6,
            x_axis_label="x",
            y_axis_label="bf(x)",
            num_graph_anchor_points=1000,
                    **kwargs)

    def construct(self):
        self.setup_axes()
        """
        BS = Bspline()
        x = np.arange(0,10,1)
        data = [1, 2, 2, 4, 4, 2, 4.2, 4, 5.5, 5.8]
        nr_splines = 5
        coef_, B, k = BS.fit(x, data, nr_splines=nr_splines).values()

        basis, knots = Bspline.basismatrix(X=np.linspace(0,10,50), nr_splines=5, l=3, knot_type="e").values()
        
        b1 = [self.get_graph(lambda x: BS.basisfunction(X=x, knots=knots, j=i+3, l=1),
                                                        x_min=0, x_max=10) for i in range(nr_splines-2)]
        b2 = [self.get_graph(lambda x: BS.basisfunction(X=x, knots=knots, j=i+3, l=2), 
                                                        x_min=0, x_max=10) for i in range(nr_splines-1)]
        b3 = [self.get_graph(lambda x: BS.basisfunction(X=x, knots=knots, j=i+3, l=3), 
                                                        x_min=0, x_max=10) for i in range(nr_splines)]
        """
        text01 = TextMobject("Basic B-splines")
        text02 = TextMobject("Of different orders")
        self.play(Write(text01))
        self.wait()
        self.play(Transform(text01, text02))
        self.wait()
        
        text_order1 = Tex("Order = 1")
        #text_order2 = Tex("Order = 2")
        #text_order3 = Tex("Order = 3")

        self.play(Transform(text01, text_order1))
        self.wait()
        """
        for bf1 in b1:
            self.add(bf1)
            self.wait()

        self.play(Transform(text01, text_order2))

        for idx, bf2 in enumerate(b2):
            if idx < len(b1):
                self.play(Transform(b1[idx], bf2))
                self.add(bf2)
                self.remove(b1[idx])
                self.wait()
            else:
                self.add(bf2)
    
        self.play(Transform(text01, text_order3))

        for idx, bf3 in enumerate(b3):
            if idx < len(b2):
                self.play(Transform(b2[idx], bf3))
                self.add(bf3)
                self.remove(b2[idx])
                self.wait()
            else:
                self.add(bf3)
                self.wait()

        data = [1, 2, 2, 4, 4, 2, 4.2, 4, 5.5, 5.8]
        for time, dat in enumerate(data):
            dot = Dot().move_to(self.coords_to_point(time, dat))
            self.add(dot)
            self.wait()
        
        """
        #l1 = Line(start=[0,0,0], end=[1, 0, 0])
        #self.play(MoveAlongPath(curve1, l1), rate_func=linear)


class SinAndCosFunctionPlot(GraphScene):
    def __init__(self, **kwargs):
        GraphScene.__init__(
            self,
            x_min=-10,
            x_max=10.3,
            num_graph_anchor_points=100,
            y_min=-1.5,
            y_max=1.5,
            graph_origin=ORIGIN,
            axes_color=GREEN,
            x_labeled_nums=range(-10, 12, 2),
            **kwargs
        )
        self.function_color = RED

    def construct(self):
        self.setup_axes(animate=False)
        func_graph = self.get_graph(np.cos, self.function_color)
        func_graph2 = self.get_graph(np.sin)
        vert_line = self.get_vertical_line_to_graph(TAU, func_graph, color=YELLOW)
        graph_lab = self.get_graph_label(func_graph, label="\\cos(x)")
        graph_lab2 = self.get_graph_label(func_graph2, label="\\sin(x)",
                            x_val=-10, direction=UP / 2)
        two_pi = MathTex(r"x = 2 \pi")
        label_coord = self.input_to_graph_point(TAU, func_graph)
        two_pi.next_to(label_coord, RIGHT + UP)
        self.add(func_graph, func_graph2, vert_line, graph_lab, graph_lab2, two_pi)