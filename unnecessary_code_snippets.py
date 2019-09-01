#  Developed by Lukas Sauer at the Heidelberg Institute of Theoretical Studies on 17.05.19, 10:01.
#  Contact: lukas.sauer@h-its.org.
#  Last modified on 17.05.19, 10:01.
#  (C) 2019. All rights reserved.


#FlagComplex.py
def get_quadratic_form_from_flags(self, flags):
    v_0 = self.ps[flags[0]]
    v_2 = self.ps[flags[1]]
    v_0_2 = line_intersection(self.ls[flags[0]], self.ls[flags[1]])


    basech = np.array([v_0, v_0_2, v_2]).T

    quadratic_form = np.array([[0, 0, 1],
                               [0, 1, 0],
                               [1, 0, 0]])
    quadratic_form = np.matmul(basech, np.matmul(quadratic_form, la.inv(basech)))
    # symmetrize matrix
    quadratic_form = 0.5*(quadratic_form+quadratic_form.T)

    self.quadratic_form = quadratic_form
    return quadratic_form
#in FlagComplex.visualize
"""
Need argument ellipse_quad=True
# Draw ellipse through point specified in quad
if ellipse_quad is not None:
    points = [self.drawps[i] for i in ellipse_quad]
    # The center is the intersection of the two lines connecting p_0 and p_2 resp. p_1 and p_3.
    center = line_intersection(connecting_line(self.ps[ellipse_quad[0]], self.ps[ellipse_quad[2]]),
                               connecting_line(self.ps[ellipse_quad[1]], self.ps[ellipse_quad[3]]))
    center = self.get_two_dimensional_point(center)
    ellipse = get_ellipse_from_points_and_center(points, center, resolution=32)
    draw_polygon(d, ellipse, col="grey")
"""
"""
if with_ellipse:
    if ellipse_flags is not None:
        quadratic_form = self.get_quadratic_form_from_flags(ellipse_flags)
    else:
        quadratic_form = self.quadratic_form
    # whatis = np.matmul(self.ps[0].T, np.matmul(quadratic_form, self.ps[0]))
    # inter = line_intersection(self.ls[0], self.ls[2])
    # whatis2 = np.matmul(self.ps[2].T, np.matmul(quadratic_form, self.ps[2]))
    # whatis_inter = np.matmul(inter.T, np.matmul(quadratic_form, inter))
    # whatis3 = np.matmul(self.ps[1].T, np.matmul(quadratic_form, self.ps[1]))
    # We see: The points p_0 and p_2 lie on the line, but the intersection point of l_0 and l_2 doesn't.
    rotation_matrix = rotate_vectors(self.projection_plane, np.array([0, 0, 1]))
    rotated_form = np.matmul(rotation_matrix, np.matmul(quadratic_form, rotation_matrix.T))
    # p0 = np.array([self.drawps[0][0], self.drawps[0][1], 1])
    # p2 = np.array([self.drawps[2][0], self.drawps[2][1], 1])
    # p3 = np.array([self.drawps[3][0], self.drawps[3][1], 1])
    # whatis_draw = np.matmul(p0.T, np.matmul(rotated_form, p0))
    # whatis_draw2 = np.matmul(p2.T, np.matmul(rotated_form, p2))
    # whatis_draw3 = np.matmul(p3.T, np.matmul(rotated_form, p3))
    ellipse = get_ellipse_from_quadratic_form(rotated_form, resolution=32)
    draw_polygon(d, [ellipse], col=["green"], wop=False)
"""

#From EuklGeometryUtility.py
def get_ellipse_from_points_and_center(points, center, resolution):
    """
    Returns a two-dimensional ellipse going through specified four points with specified center.
    Resolution is the number of points on the ellipse (times four) that will be returned.

    WARNING! Until now, I didn't get this to work...

    :param points: a list of four two-dimensional numpy arrays
    :param center: a two-dimensional numpy array
    :return: a list of two-dimensional numpy arrays (points on the ellipse, evenly spanning 360 degrees)
    """
    # The most general form of an ellipse is the normal form of the quadric, i.e.
    #
    # Ax^2 + Bxy + Cy^2 + Dx + Ey = 1 [Equation I]
    #
    # As we have four points, we can solve this linear equation depending on the last variable E.
    # We do this by this system as matrix M  with rows (x^ 2, xy, y^2, x) for each point in points. Then we have
    #
    # v = (A, B, C, D)^T = - M^-1 * E * (y_i)_i + M^-1 * (1)_i  = E * c + d
    #
    # where (y_i)_i is the vector of y-components of the points p_i, and (1)_i is a vector consisting only of ones.
    # and c= M^-1*(y_i)_i and d = M^-1 * (1)_i.
    #
    # For determining E, we use the fact that the variables of the center satisfy:
    #
    # x_center = (2CD - BE)/(B^2 - 4AC) and y_center = (2AE - BD)/(B^2 - 4AC) [Equations II and III]
    #
    # Both of them are quadratic equations in E, so their common solution is the variable E.

    M = np.zeros((4, 4))
    c = np.zeros((4, 1))
    d = np.zeros((4, 1))

    for i in range(4):
        M[i] = np.array([points[i][0]**2, points[i][0]*points[i][1], points[i][1]**2, points[i][0]])
        c[i] = points[i][1]
        d[i] = 1

    c = np.matmul(la.inv(M), c)
    d = np.matmul(la.inv(M), d)

    # The quadratic equations II and III coeff[0]*E**2+coeff[1]*E +coeff[2] have the coefficients below,
    # and can be solved with the built-in python module.

    coeff_x = np.array([0, 0, 0])
    coeff_x[0] = center[0]*(c[1]**2 - 4*c[0]*c[2]) - 2*c[2]*c[3] + c[1]
    coeff_x[1] = center[0]*(2*c[1]*d[1] - 4*(c[0]*d[2] + c[2]*d[0])) - 2*(c[2]*d[3] + c[3]*d[2]) + d[1]
    coeff_x[2] = center[0]*(d[1]**2 - 4*d[0]*d[2]) - 2*d[2]*d[3]

    coeff_y = np.array([0, 0, 0])
    coeff_y[0] = center[1]*(c[1]**2 - 4*c[0]*c[2]) - 2*c[0] + c[1]*c[3]
    coeff_y[1] = center[1]*(2*c[1]*d[1] - 4*(c[0]*d[2]+c[2]*d[0])) - 2*d[0] + c[1]*d[3] + c[3]*d[1]
    coeff_y[2] = center[1]*(d[1]**2 - 4*d[0]*d[2]) + d[1]*d[3]

    x_roots = np.roots(coeff_x)
    y_roots = np.roots(coeff_y)

    common = set(x_roots).intersection(y_roots)

    assert len(common) == 1, "Warning! The number of common roots is not 1, but " + str(len(common)) + "!"

    E = common.pop()

    # The vector v=(A,B,C,D)^T

    v = E*c + d
    A = v[0]
    B = v[1]
    C = v[2]
    D = v[3]
    # The canonical form is also written as Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0, so in our case:
    F = -1

    # Now that we have calculated the coefficients (A, B, C, D, E) of the normal form, we can use them to
    # calculate semi-major and semi-minor axis and the angle that the ellipse is rotated to the canonical form,
    # i.e. the angle between the semi-major axis and the horizontal axis of the coordinate system.
    #
    # Check the article https://en.wikipedia.org/wiki/Ellipse for understanding the formulas.
    # The semi_major and semi_minor are called a resp. b in the article.

    # The two axes differ only by one minus sign.                        ---> over here
    semi_major = - np.sqrt(2 * (A*E**2 + C*D**2 - B*D*E + (B**2 - 4*A*C)*F)(A + C + np.sqrt((A - C)**2 + B**2))) / (B**2 - 4*A*C)
    semi_minor = - np.sqrt(2 * (A*E**2 + C*D**2 - B*D*E + (B**2 - 4*A*C)*F)(A + C - np.sqrt((A - C)**2 + B**2))) / (B**2 - 4*A*C)

    if B == 0:
        if A < C:
            angle = 0
        elif C > A:
            # 90 degrees
            angle = np.pi/2
        else:
            Exception("Warning! A equals C, so the ellipse is degenerate!")
    else:
        s = (C - A - np.sqrt((A - C)**2 + B**2))/B
        angle = np.arctan(s)

    # Now we are ready to generate the ellipse. We generate an ellipse with center (0,0) and semi_major and semi_minor
    # oriented along the coordinate axes. Afterwards, we shift it to our center, and rotate it with the respective angle.

    ellipse = []

    for i in range(4*resolution):
        phi = np.pi/(2*resolution)*i
        # Ellipse coordinates
        point = np.array([semi_major*np.cos(phi), semi_minor*np.sin(phi)])
        # Rotation
        R = rotate_by_angle_2dim(angle)
        point = R*point
        # Translation
        point = point + center
        ellipse.append(point)

    return ellipse

def get_ellipse_from_quadratic_form(matrix, resolution ,hyperbola_length = 5):

    a = matrix[0][0]
    b = matrix[0][1]*2
    c = matrix[1][1]
    d = matrix[2][0]*2
    e = matrix[1][2]*2
    f = matrix[2][2]

    center = np.array([(2 * c * d - b * e) / (b ** 2 - 4 * a * c), (2 * a * e - b * d) / (b ** 2 - 4 * a * c)])

    #print("Ellipse imaginary if 0< " + str(la.det(matrix)*(a+c)))
    # assert la.det(matrix)*(a+c) < 0, "The ellipse is imaginary or a point!"

    # return get_ellipse_from_conic_section_equation([a, b, c, d, e, f], resolution)


    # Get the eigenvalues of the matrix A_{33}
    eigenvalues, eigenvectors = la.eig(matrix[:2, :2])
    det_A_Q = la.det(matrix)
    det_A_33 = la.det(matrix[:2, :2])

    is_ellipse = np.sign(eigenvalues[0]) == np.sign(eigenvalues[1])

    axes_squared = np.array([1/x * (-1) * det_A_Q / det_A_33 for x in eigenvalues])
    idx = axes_squared.argsort()[::-1]
    axes_squared = axes_squared[idx]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    conic = []

    if is_ellipse:
        conic = get_ellipse_from_axes(np.sqrt(axes_squared[0]), np.sqrt(axes_squared[1]), resolution)
    else: # is hyperbola
        conic = get_hyperbola_from_axes(np.sqrt(axes_squared[0]), np.sqrt(-axes_squared[1]), resolution, hyperbola_length)

    conic = [np.matmul(eigenvectors, point) + center for point in conic]

    return conic

def get_ellipse_from_conic_section_equation(coeff, resolution):
    [a, b, c, d, e, f] = coeff

    assert a*c-b**2/4 < 0,\
        "Attention! The coefficients a to f define a hyperbola or a parabola, but not an ellipse."

    # Coefficient normalizing factor
    #q = 64*(f*(4*a*c-b**2) - a*e**2 + b*d*e - c*d**2)/(4*a*c - b**2)**2
    # Distance between center and focal point
    #s = 0.25 * np.sqrt(np.abs(q)*np.sqrt(b**2+(a-c)**2))
    # Semi-major axis
    #semi_major = 0.125 * np.sqrt(2*np.abs(q)*np.sqrt(b**2+(a-c)**2)-2*q*(a+c))
    # Semi-minor axis
    #semi_minor = np.sqrt(semi_major**2 + s**2)
    # Now that we have calculated the coefficients (a, b, c, d, e) of the normal form, we can use them to
    # calculate semi-major and semi-minor axis and the angle that the ellipse is rotated to the canonical form,
    # i.e. the angle between the semi-major axis and the horizontal axis of the coordinate system.
    #
    # Check the article https://en.wikipedia.org/wiki/Ellipse for understanding the formulas.
    # The semi_major and semi_minor are called a resp. b in the article.

    # The two axes differ only by one minus sign.                        ---> over here
    semi_major = - np.sqrt(2 * (a * e ** 2 + c * d ** 2 - b * d * e + (b ** 2 - 4 * a * c) * f)*(
        a + c + np.sqrt((a - c) ** 2 + b ** 2))) / (b ** 2 - 4 * a * c)
    step1 = 2 * (a * e ** 2 + c * d ** 2 - b * d * e + (b ** 2 - 4 * a * c) * f)
    step2 = (a + c - np.sqrt((a - c) ** 2 + b ** 2))
    step3 = np.sqrt(step1*step2)
    semi_minor = - step3 / (b ** 2 - 4 * a * c)

    center = np.array([(2*c*d - b*e)/(b**2 - 4*a*c), (2*a*e - b*d)/(b**2 - 4*a*c)])
    """
    if q*a - q*c == 0 and q*b == 0:
        angle = 0
    if q*a - q*c == 0 and q*b > 0:
        angle = 0.25*np.pi
    if q*a - q*c == 0 and q*b < 0:
        angle = 0.75*np.pi
    if q*a - q*c > 0 and q*b >= 0:
        angle = 0.5*np.arctan(b*a - c)
    if q*a - q*c > 0 and q*b < 0:
        angle = 0.5*np.arctan(b*a - c)+np.pi
    if q*a - q*c < 0:
        angle = 0.5*np.arctan(b*a - c)+0.5*np.pi
    """


    if b == 0:
        if a < c:
            angle = 0
        elif c > a:
            # 90 degrees
            angle = np.pi / 2
        else:
            Exception("Warning! a equals c, so the ellipse is degenerate!")
    else:
        s = (c - a - np.sqrt((a - c) ** 2 + b ** 2)) / b
        angle = np.arctan(s)

    return get_ellipse_from_axes(semi_major,semi_minor, center, angle, resolution)



def get_ellipse_from_axes(semi_major, semi_minor, center=np.array([0,0]), angle=0, resolution=32):
    # # Now we are ready to generate the ellipse. We generate an ellipse with center (0,0) and semi_major and semi_minor
    # # oriented along the coordinate axes. Afterwards, we shift it to our center, and rotate it with the respective angle.

    ellipse = []

    for i in range(4 * resolution):
        phi = np.pi / (2 * resolution) * i
        # Ellipse coordinates
        point = np.array([semi_major * np.cos(phi), semi_minor * np.sin(phi)])
        # Rotation
        if angle != 0:
            R = rotate_by_angle_2dim(angle)
            point = np.matmul(R,  point)
        # Translation
        point = point + center
        ellipse.append(point)

    return ellipse

def get_hyperbola_from_axes(semi_major, semi_minor, resolution=32, hyperbola_length=5):
    # # Now we are ready to generate the hyperbola. We generate an hyperbola with center (0,0) and semi_major and semi_minor
    # # oriented along the coordinate axes. Afterwards, we shift it to our center, and rotate it with the respective angle.

    hyperbola = []

    for i in range(2 * resolution):
        t = (-  0.5 + i / (2 * resolution))*hyperbola_length
        # Hyperbolic coordinates
        point = np.array([semi_major * np.cosh(t), semi_minor * np.sinh(t)])
        hyperbola.append(point)

    for i in range(2*resolution):
        t = (0.5 - i / (2 * resolution)) * hyperbola_length
        # Hyperbolic coordinates
        point = np.array([-semi_major * np.cosh(t), semi_minor * np.sinh(t)])
        hyperbola.append(point)

    return hyperbola








