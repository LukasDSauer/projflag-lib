#  Developed by Lukas Sauer at the Heidelberg Institute of Theoretical Studies on 3/29/19, 10:37 AM.
#  Contact: lukas.sauer@h-its.org.
#  Last modified on 3/29/19, 10:37 AM.
#  (C) 2019. All rights reserved.


class Model3DMaker:

    def __init__(self, transformation_width, get_shape, transform_shape, transform_args, steps=100, dist_per_step=0.06, filename="model3d.obj"):
        self.steps = steps
        self.transformation_width = transformation_width
        self.filename = filename
        self.get_shape = get_shape
        self.transform_shape = transform_shape
        self.transform_args = transform_args
        self.step_width = self.transformation_width/self.steps
        self.dist_per_step = dist_per_step
        self.n = len(self.get_shape())

    def make_model(self, filename=""):
        if filename == "":
            file = open(self.filename, "w")
        else:
            file = open(filename, "w")

        for r in range(self.steps):
            args = self.transform_args
            args["t"] = self.step_width
            self.transform_shape(**args)
            points = self.get_shape()

            for p in points:
                file.write("v " + str(p[0]) + " " + str(p[1]) + " " + str(r * self.dist_per_step) + "\n")

        for r in range(self.steps - 1):
            for i in range(self.n):
                file.write("f " + str(i + 1 + r * 3) + " " + str(i + 2 + r * 3) + " " + str(i + 1 + self.n + r * 3) + "\n")
                file.write("f " + str(i + 2 + r * 3) + " " + str(i + 2 + self.n + r * 3) + " " + str(i + 1 + self.n + r * 3) + "\n")

        file.close()










