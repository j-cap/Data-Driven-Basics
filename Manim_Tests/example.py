
from manim import *

# to run: 
# manim -p -ql example.py SquareToCircle 

class SquareToCircle(Scene):
    def construct(self):
        circle = Circle()
        square = Square()
        square.flip(RIGHT)
        square.rotate(-3 * TAU / 8)
        circle.set_fill(PINK, opacity=0.5)

        self.play(Create(square))
        self.play(Transform(square, circle))
        self.play(FadeOut(square))
        self.play(Create(circle))
        self.play(Transform(circle, square))
        self.play(FadeOut(circle))