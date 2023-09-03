package upb.problem;

import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.problem.AnalyticalProblem;

import java.util.Arrays;

import static java.lang.Math.*;
import static upb.problem.FunctiiBobine.campHaxa_bobina_supSInfSub;
import static upb.problem.FunctiiBobine.campHoriunde_bobina_supSpElem;
import static upb.problem.MatlabFunctions.linspace;

public class LoneyMO extends AbstractProblem implements AnalyticalProblem {

    private double[] HBobinaFixa_0;
    private double h = 12e-2;
    private double z0 = 2.5e-3;
    private double zs = 24e-2;
    private double r1 = 1.1e-2;
    private double r2 = 2.9e-2;
    private int NOSOL = 4;
    private double r3 = 3e-2;
    private double r4 = 3.6e-2;
    private int NOSP = 1000;
    private int NOPAXA = 6;
    private int ni = 1;
    private double zmin = -h / 2;
    private double[] vectorz;
    private double CrtTotBobiFix = ni * h * NOSOL;
    private double tol = ulp(1);
    private double rs = 10.8e-2;
    private double[][] coord = new double[][]{new double[]{rs, 0}};
    private double[][] HBobinaFixa_1;

    public LoneyMO() {
        super(2, 2);
        vectorz = Arrays.copyOf(linspace(0, z0, NOPAXA), NOPAXA + 1);
        vectorz[NOPAXA] = zs;
        HBobinaFixa_0 = campHaxa_bobina_supSInfSub(r1, r2, zmin, h, CrtTotBobiFix, vectorz, NOSOL);
        HBobinaFixa_1 = campHoriunde_bobina_supSpElem(r1, r2, zmin, h, CrtTotBobiFix, coord, NOSOL, NOSP, tol);
    }

    @Override
    public void evaluate(Solution solution) {
        double[] x = EncodingUtils.getReal(solution);
        double s = x[0];
        double l = x[1];
        if (x[0] < 0 || x[0] > 1.2 || x[1] < 0 || x[1] > 0.7) throw new RuntimeException();
        double CrtTotBobiCor = ni * l * NOSOL;
        double zmin = s / 2;
        double zmin2 = -s / 2 - l;

        double[] HBobinaCor_0 = campHaxa_bobina_supSInfSub(r3, r4, zmin, zmin2, l, CrtTotBobiCor, vectorz, NOSOL);
        double[][] HbobinaCor_1 = campHoriunde_bobina_supSpElem(r3, r4, zmin, zmin2, l, CrtTotBobiCor, coord, NOSOL, NOSP, tol);

        double Hmin = 0;
        double Hmax = 0;
        for (int i = 0; i < vectorz.length; i++) {
            HBobinaCor_0[i] += HBobinaFixa_0[i];
            if (i == 0) {
                Hmin = HBobinaCor_0[0];
                Hmax = HBobinaCor_0[0];
                continue;
            }
            if (i == vectorz.length - 1) {
                continue;
            }
            if (HBobinaCor_0[i] > Hmax) {
                Hmax = HBobinaCor_0[i];
            }
            if (HBobinaCor_0[i] < Hmin) {
                Hmin = HBobinaCor_0[i];
            }

        }
        for (int i = 0; i < coord.length; i++) {
            HbobinaCor_1[i][0] += HBobinaFixa_1[i][0];
            HbobinaCor_1[i][1] += HBobinaFixa_1[i][1];
        }

        double H0 = HBobinaCor_0[0];
        double D = (Hmax - Hmin) / H0;

        double H1 = abs(HBobinaCor_0[HBobinaCor_0.length - 1]);
        double H2 = abs(HbobinaCor_1[0][1]);

        double Delta = max(H1, H2) / H0;
        solution.setObjectives(new double[]{D, Delta});
    }

    @Override
    public Solution newSolution() {
        Solution newSolution = new Solution(2, 2);
        newSolution.setVariable(0, new RealVariable(0, 0, 1.2));
        newSolution.setVariable(1, new RealVariable(0, 0, 0.7));
        return newSolution;
    }

    @Override
    public Solution generate() {
        Solution solution = newSolution();
        double s = PRNG.nextDouble(0.0, 1.2);
        double l = PRNG.nextDouble(0.0, 0.7);

        EncodingUtils.setReal(solution.getVariable(0), s);
        EncodingUtils.setReal(solution.getVariable(1), l);

        evaluate(solution);
        return solution;
    }

    public double getH() {
        return h;
    }

    public void setH(double h) {
        this.h = h;
    }

    public double getZ0() {
        return z0;
    }

    public void setZ0(double z0) {
        this.z0 = z0;
    }

    public double getZs() {
        return zs;
    }

    public void setZs(double zs) {
        this.zs = zs;
    }

    public double getR1() {
        return r1;
    }

    public void setR1(double r1) {
        this.r1 = r1;
    }

    public double getR2() {
        return r2;
    }

    public void setR2(double r2) {
        this.r2 = r2;
    }

    public int getNOSOL() {
        return NOSOL;
    }

    public void setNOSOL(int NOSOL) {
        this.NOSOL = NOSOL;
    }

    public double getR3() {
        return r3;
    }

    public void setR3(double r3) {
        this.r3 = r3;
    }

    public double getR4() {
        return r4;
    }

    public void setR4(double r4) {
        this.r4 = r4;
    }

    public int getNOSP() {
        return NOSP;
    }

    public void setNOSP(int NOSP) {
        this.NOSP = NOSP;
    }

    public int getNOPAXA() {
        return NOPAXA;
    }

    public void setNOPAXA(int NOPAXA) {
        this.NOPAXA = NOPAXA;
    }

    public int getNi() {
        return ni;
    }

    public void setNi(int ni) {
        this.ni = ni;
    }

    public double getZmin() {
        return zmin;
    }

    public void setZmin(double zmin) {
        this.zmin = zmin;
    }

    public double[] getVectorz() {
        return vectorz;
    }

    public void setVectorz(double[] vectorz) {
        this.vectorz = vectorz;
    }

    public double getCrtTotBobiFix() {
        return CrtTotBobiFix;
    }

    public void setCrtTotBobiFix(double crtTotBobiFix) {
        CrtTotBobiFix = crtTotBobiFix;
    }

    public double[] getHBobinaFixa_0() {
        return HBobinaFixa_0;
    }

    public void setHBobinaFixa_0(double[] HBobinaFixa_0) {
        this.HBobinaFixa_0 = HBobinaFixa_0;
    }

    public double getTol() {
        return tol;
    }

    public void setTol(double tol) {
        this.tol = tol;
    }

    public double getRs() {
        return rs;
    }

    public void setRs(double rs) {
        this.rs = rs;
    }

    public double[][] getCoord() {
        return coord;
    }

    public void setCoord(double[][] coord) {
        this.coord = coord;
    }

    public double[][] getHBobinaFixa_1() {
        return HBobinaFixa_1;
    }

    public void setHBobinaFixa_1(double[][] HBobinaFixa_1) {
        this.HBobinaFixa_1 = HBobinaFixa_1;
    }
}
