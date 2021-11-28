from math import cos, sin, tan
from matplotlib.pyplot import plot
import numpy as np

class Contact(object):
    """
       Rectangle Contact built from a center point and 
    """
    def __init__(self, center, yaw=0.0, half_length=0.6, half_width=0.4) -> None:
        super().__init__()
        self.center = center
        self.half_length = half_length
        self.half_width = half_width
        self.yaw = yaw

        self.low_bound = np.array([[self.center[0]+tan(self.yaw)*self.center[1]-self.half_length/abs(cos(self.yaw))],
                                   [self.center[1]-tan(self.yaw)*self.center[0]-self.half_width/abs(cos(self.yaw))]])
        self.constraint = np.array([[1, tan(self.yaw), 0, 0, 0, 0],
                                    [-tan(self.yaw), 1, 0, 0, 0, 0]])
        self.up_bound = np.array([[self.center[0]+tan(self.yaw)*self.center[1]+self.half_length/abs(cos(self.yaw))],
                                  [self.center[1]-tan(self.yaw)*self.center[0]+self.half_width/abs(cos(self.yaw))]])
        self.rotation_matrix = np.array([[cos(self.yaw), sin(self.yaw)],
                                         [-sin(self.yaw), cos(self.yaw)]])
        self.visual_matrix = self.center + np.array([[self.half_length, self.half_width],
                                                        [self.half_length,-self.half_width],
                                                        [-self.half_length,-self.half_width],
                                                        [-self.half_length,self.half_width],
                                                        [self.half_length, self.half_width]]) @ self.rotation_matrix

    def inside(self, state):
        if (self.constraint @ state >= self.low_bound).all() and (self.constraint @ state <= self.up_bound).all():
            print("in")
        else:
            print("out")
        
        self.show(state)

    def show(self, p):
        import matplotlib.pyplot as plt
        plt.plot(self.visual_matrix[:, 0], self.visual_matrix[:, 1])
        plt.scatter(p[0], p[1])
        plt.axis("equal")
        plt.show()

if __name__ == "__main__":
    center = np.array([0, 0])
    yaw = np.math.pi / 4.
    a = Contact(center, yaw)
    state = np.random.randn(6)

    a.inside(state)
