from manim import *

from stareg.bspline import Bspline
import numpy as np


class PBS(GraphScene):
    def __init__(self, **kwargs):
        GraphScene.__init__(
            self, 
            y_min=-1.2,
            y_max=1.2,
            x_min=0.,
            x_max=1., 
            **kwargs
        )

    def construct(self):
        self.setup_axes(animate=True)

        BS = Bspline()
        X = np.random.uniform(0,1,10)
        Y = np.sin(2*np.pi*X) + np.random.normal(scale=0.05, size=10)
       
        self.setup_axes()
        dot_collection = VGroup()
        for time, val in enumerate(Y):
            dot = Dot().move_to(self.coords_to_point(X[time], val))
            self.add(dot)
            dot_collection.add(dot)


        B, knots = BS.basismatrix(np.linspace(0,1,1000)).values()
        

        graph = self.get_graph(
            lambda x: float(BS.basisfunction(np.asarray(x), knots, 5, 1)), x_min=-0.2, x_max=1.2)
        graph2 = self.get_graph(
            lambda x: float(BS.basisfunction(np.asarray(x), knots, 6, 2)), x_min=-0.2, x_max=1.2)
        graph3 = self.get_graph(
            lambda x: float(BS.basisfunction(np.asarray(x), knots, 7, 3)), x_min=-0.2, x_max=1.2)
        
        graphs_linear = []
        for i in range(0,B.shape[1]+1):
            graphs_linear.append(self.get_graph(lambda x: float(
                BS.basisfunction(np.asarray(x), knots, i, 1)), x_min=0, x_max=1))

        #self.play(Create(graph))      # animate the creation of the square

        for g in graphs_linear:
            self.play(Create(g))
        #path = self.get_graph(lambda x: 0, x_min=0, x_max=0.5)
        #self.play(ShowCreation(path))
        #self.play(MoveAlongPath(graph_copy, path))
        #self.play(Transform(graph, graph2), runtime=2) # interpolate the square into the circle
        #self.play(Transform(graph2, graph3), runtime=2) # interpolate the square into the circle

        


class Test(GraphScene):
    def __init__(self, **kwargs):
        GraphScene.__init__(
            self, 
            #y_min=-2,
            #y_max=25,
            #x_min=-3,
            #x_max=10, 
            **kwargs
        )

    def get_points_from_coords(self, coords):
        return [self.coords_to_point(px, py) for px, py in coords]
    def get_dots_from_coords(self, coords, radius=0.1):
        points = self.get_points_from_coords(coords)
        dots = VGroup(*[Dot(radius=radius).move_to([px,py, pz]) for px, py, pz in points])
        return dots

    def construct(self):
        self.setup_axes(animate=True)
        #path = self.get_graph(lambda x: (x-2)**2, x_min=2, x_max=-3)
        #location = self.coords_to_point(2,0) #location: Point
        #dot = Dot(location) #dot: Dot
        #self.play(ShowCreation(path), ShowCreation(dot))
        #self.play(MoveAlongPath(dot, path))

        x = [0,1,2,3,4,5,6,7]
        y = [0,1,4,9,16,25,20,10]
        coords = [[px, py] for px, py in zip(x,y)]
        #points = self.get_points_from_coords(coords)
        dots = self.get_dots_from_coords(coords, 0.1)
        self.add(dots)



class PlotBSpline(Scene):
    # contributed by heejin_park, https://infograph.tistory.com/230
    def construct(self):
        self.show_axis()
        #self.show_circle()
        #self.move_dot_and_draw_curve()
        #self.wait()

    def show_axis(self):
        x_start = np.array([-1,0,0])
        x_end = np.array([6,0,0])

        y_start = np.array([0,-2,0])
        y_end = np.array([0,2,0])

        x_axis = Line(x_start, x_end)
        y_axis = Line(y_start, y_end)

        self.add(x_axis, y_axis)

        self.origin_point = np.array([0,0,0])
        self.curve_start = np.array([0,0,0])

    def show_circle(self):
        circle = Circle(radius=1)
        circle.move_to(self.origin_point)
        self.add(circle)
        self.circle = circle

    def move_dot_and_draw_curve(self):
        orbit = self.circle
        origin_point = self.origin_point

        dot = Dot(radius=0.08, color=YELLOW)
        dot.move_to(orbit.point_from_proportion(0))
        self.t_offset = 0
        rate = 0.25

        def go_around_circle(mob, dt):
            self.t_offset += (dt * rate)
            # print(self.t_offset)
            mob.move_to(orbit.point_from_proportion(self.t_offset % 1))

        def get_line_to_circle():
            return Line(origin_point, dot.get_center(), color=BLUE)

        def get_line_to_curve():
            x = self.curve_start[0] + self.t_offset * 4
            y = dot.get_center()[1]
            return Line(dot.get_center(), np.array([x,y,0]), color=YELLOW_A, stroke_width=2 )


        self.curve = VGroup()
        self.curve.add(Line(self.curve_start,self.curve_start))
        def get_curve():
            last_line = self.curve[-1]
            x = self.curve_start[0] + self.t_offset * 4
            y = dot.get_center()[1]
            new_line = Line(last_line.get_end(),np.array([x,y,0]), color=YELLOW_D)
            self.curve.add(new_line)

            return self.curve

        dot.add_updater(go_around_circle)

        origin_to_circle_line = always_redraw(get_line_to_circle)
        dot_to_curve_line = always_redraw(get_line_to_curve)
        sine_curve_line = always_redraw(get_curve)

        self.add(dot)
        self.add(orbit, origin_to_circle_line, dot_to_curve_line, sine_curve_line)
        self.wait(8.5)

        dot.remove_updater(go_around_circle)