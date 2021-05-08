
from manim import *

class SinAndCosFunctionPlot(GraphScene):
    arguments = {
        "x_min": -10,
        "x_max": 10.3,
        "num_graph_anchor_points": 100,
        "y_min": -1.5,
        "y_max": 1.5,
        "graph_origin": ORIGIN,
        "function_color": RED,
        "axes_color": GREEN,
        "x_labeled_nums": range(-10, 12, 2),
    }

    def construct(self):
        self.setup_axes(animate=False)
        func_graph = self.get_graph(np.cos, self.arguments["function_color"])
        func_graph2 = self.get_graph(np.sin)
        vert_line = self.get_vertical_line_to_graph(TAU, func_graph, color=YELLOW)
        graph_lab = self.get_graph_label(func_graph, label="\\cos(x)")
        graph_lab2 = self.get_graph_label(func_graph2, label="\\sin(x)",
                            x_val=-10, direction=UP / 2)
        two_pi = MathTex(r"x = 2 \pi")
        label_coord = self.input_to_graph_point(TAU, func_graph)
        two_pi.next_to(label_coord, RIGHT + UP)
        self.add(func_graph, func_graph2, vert_line, graph_lab, graph_lab2, two_pi)