#  Developed by Lukas Sauer at the Heidelberg Institute of Theoretical Studies on 10/30/19, 5:20 PM.
#  Contact: lukas.sauer@h-its.org.
#  Last modified on 10/30/19, 4:41 PM.
#  (C) 2019. All rights reserved.

#  Developed by Lukas Sauer at the Heidelberg Institute of Theoretical Studies on 4/29/19, 11:38 AM.
#  Last modified on 4/29/19, 11:38 AM.
#  (C) 2019. All rights reserved.

import numpy as np

from scipy.spatial import ConvexHull


import string

from .ProjGeometryUtility import get_reflection, connecting_line



class FlagTesselator:
    def __init__(self, flagcomplex, steps=3):
        """

        :param flagcomplex: an object of the flagcomplex class
        :param steps: an integer, the number of steps that tesselation is iterated
        """
        # The flagcomplex object containing all the flags
        self.flagcomplex = flagcomplex
        # The number of flags
        self.n = self.flagcomplex.n
        # The number of iterations of the tesselation
        self.steps = steps
        # The initial polygon
        self.initial = None
        # This will be the connecting lines between p_i-1 and p_i. Along these lines, the reflections will be performed.
        self.connecting_ls = None
        # The intersection points between the flags
        self.qs = None
        # The reflection maps
        self.reflections = None
        # The different reflections generate a group. This group can be represented by words, where each reflection
        # is represented by a different letter.
        self.words = None
        self.letters = None
        # The important maps for the tesselation are represented by words up to a certain length, the number of steps.
        self.maps = None

    def refresh_setup(self):
        """
        Makes sure all the coordinates in the flag complex are up to date. Then it computes the necessary reflections and
        the words up to the length = self.steps. Finally, it computes the maps that are represented by those words.

        :return:
        """
        self.flagcomplex.refresh_setup()
        self.initial = np.array([p for p in self.flagcomplex.ps])
        self.qs = self.flagcomplex.qs
        # The lines we want to reflect along are the connecting lines of points p_i and p_i-1.
        self.connecting_ls = [connecting_line(self.flagcomplex.ps[i - 1], self.flagcomplex.ps[i]) for i in range(self.n)]
        self.compute_reflections()
        self.compute_words()
        self.compute_maps_from_words()


    def generate_tesselation(self):
        """
        This function performs the tesselation. It returns a two-dimensional projection of the initial polygon, the
        convex hull of the tesselation, and finally a dictionary of tiles (the reflections of the initial polygon,
        indexed by the words which generate them)

        :return: a triple of objects, the first is a list of numpy arrays, the second as well, the third is a dictionary
        of lists of numpy arrays
        """
        self.refresh_setup()
        tiles = dict()
        # Map the initial polygon according to the different words.
        for word in self.words:
            tiles[word] = [self.map_point_to_two_dim(self.maps[word], y) for y in self.initial]
        # Compute the convex hull.
        maxlwords = np.array([tiles[x] for x in self.words if len(x) == self.steps])
        maxlwords = maxlwords.reshape(-1, 2)
        cv = ConvexHull([x[:2] for x in maxlwords])
        hull_points = cv.vertices
        hull = [maxlwords[h] for h in hull_points]

        return [self.flagcomplex.get_two_dimensional_point(x) for x in self.initial], hull, tiles

    def map_point_to_two_dim(self, m, x):
        """
        Transforms a projective point with the projective transformation m and then projects it to two dimensions.

        :param m: a projective transformation (3x3 dim np-array)
        :param x: a projective point (three dim np-array)
        :return: a two-dim np-array
        """
        point = np.matmul(m, x)
        return self.flagcomplex.get_two_dimensional_point(point)

    def set_steps(self, s):
        """
        Set the number of steps you want to iterate.

        :param s:
        :return:
        """
        self.steps = s

    def compute_reflections(self):
        """
        Computes the reflections along all connection lines with the respective intersection point.

        :return:
        """
        self.reflections = np.array([get_reflection(self.connecting_ls[i], self.qs[i]) for i in range(self.n)])

    def compute_words(self):
        """
        Computes all the relevant words up to the length self.steps. Words with double letters are irrelevant because
        the letters cancel out.

        :return:
        """
        letters = list(string.ascii_lowercase)[:self.n]
        words = list(string.ascii_lowercase)[:self.n]
        for i in range(1, self.steps):
            tmp = words.copy()
            for w in tmp:
                for l in letters:
                    # as l*l = e, we don't need double letters in our words.
                    if w[-1] != l:
                        words.append(w + l)
        self.words = words
        self.letters = letters

    def compute_maps_from_words(self):
        """
        Computes the respective projective transformation from each word.

        :return:
        """
        self.maps = dict()

        reflection_by_letter = dict()
        for i in range(self.n):
            reflection_by_letter[self.letters[i]] = self.reflections[i]

        for word in self.words:
            current = np.eye(3)
            for w in reversed(word):
                current = np.matmul(reflection_by_letter[w], current)
            self.maps[word] = current
