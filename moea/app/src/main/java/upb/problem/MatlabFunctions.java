package upb.problem;

import static java.lang.Math.*;

class MatlabFunctions {

    static double[] linspace(double x1, double x2, int n) {
        double spacing = (x2 - x1) / (n - 1);
        double[] x = new double[n];
        x[n - 1] = x2;
        for (int i = 0; i < n - 1; i++) {
            x[i] = x1 + i * spacing;
        }
        return x;
    }

    /**
     * Computes the complete elliptic integral to accuracy tol.
     *
     * @param tol accuracy, increase its value for a less accurate but more quickly computed answer.
     * @return complete elliptic integral of the first and second kind
     */
    static double[] ellipke(double k, double tol) {
        double a1 = 1;
        double a2 = 1 - k;
        double b1 = abs(sqrt(a2));
        double b2 = abs(sqrt(a2));
        if (b1 == 0.0) return new double[]{0, 0};

        double pio2 = PI / 2;

        double
                c1 = 1, c2 = 1,
                d1 = 1, d2 = 1,
                e1 = b1, e2 = b2,
                f1 = 1, f2 = 1,
                g1 = d1, g2 = d2,
                h1, h2;

        a1 = a1 / c1;
        a2 = a2 / c2;
        d1 = a1 / c1 + d1;
        d2 = a2 / c2 + d2;
        h1 = e1 / c1;
        h2 = e2 / c2;
        a1 = 2 * (g1 * h1 + a1);
        a2 = 2 * (g2 * h2 + a2);
        c1 = c1 + h1;
        c2 = c2 + h2;
        h1 = f1;
        h2 = f2;
        f1 = f1 + b1;
        f2 = f2 + b2;

        boolean abs1 = abs(h1 - b1) > h1 * tol;
        boolean abs2 = abs(h2 - b2) > h2 * tol;

        while (abs1 || abs2) {
            if (abs1) {
                b1 = 2 * sqrt(e1);
                e1 = b1 * f1;
                g1 = d1;
                d1 = a1 / c1 + d1;
                h1 = e1 / c1;
                a1 = 2 * (g1 * h1 + a1);
                c1 = c1 + h1;
                h1 = f1;
                f1 = f1 + b1;
            }
            if (abs2) {
                b2 = 2 * sqrt(e2);
                e2 = b2 * f2;
                g2 = d2;
                d2 = a2 / c2 + d2;
                h2 = e2 / c2;
                a2 = 2 * (g2 * h2 + a2);
                c2 = c2 + h2;
                h2 = f2;
                f2 = f2 + b2;
            }
            abs1 = abs(h1 - b1) > h1 * tol;
            abs2 = abs(h2 - b2) > h2 * tol;
        }

        return new double[]{pio2 * (d1 * f1 + a1) / (f1 * (f1 + c1)), pio2 * (d2 * f2 + a2) / (f2 * (f2 + c2))};
    }

}
