package upb.algorithm;

import org.moeaframework.algorithm.pso.OMOPSO;
import org.moeaframework.analysis.sensitivity.EpsilonHelper;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;

import java.util.Random;

import static java.lang.Math.abs;

// version 6
public class OMOACORV3 extends OMOPSO {

    private static final Random rand = new Random();
    int numberOfVariables;
    int numarFurnici;
    float vitezaDeConvergenta;
    double[] deviatiileStandard;

    public OMOACORV3(Problem problem) {
        this(
                problem,
                100,
                100,
                new double[]{EpsilonHelper.getEpsilon(problem)},
                1.0 / problem.getNumberOfObjectives(),
                0.01,
                100,
                0.01f,
                100
        );
    }

    public OMOACORV3(Problem problem, int swarmSize, int leaderSize, double[] epsilons, double mutationProbability, double mutationPerturbation, int maxIterations, float vitezaDeConvergenta, int numarFurnici) {
        super(problem, swarmSize, leaderSize, epsilons, mutationProbability, mutationPerturbation, maxIterations);
        numberOfVariables = problem.getNumberOfVariables();
        this.vitezaDeConvergenta = vitezaDeConvergenta;
        deviatiileStandard = new double[swarmSize];
        this.numarFurnici = numarFurnici;
    }

    @Override
    protected void iterate() {

        constructieSolutii();

        mutate();

        evaluateAll(particles);

        updateLocalBest();

        leaders.addAll(particles);
        leaders.update();

        if (archive != null) {
            archive.addAll(particles);
        }
    }

    private void constructieSolutii() {
        for (int i = 0; i < particles.length; i++) {
            calculeazaDeviatiileStandard(i);
            particles[i] = construiesteSolutieNoua(i);
        }
    }

    private Solution construiesteSolutieNoua(int indexArhiva4Constructie) {
        Solution solutieDeStart = particles[indexArhiva4Constructie];
        Solution solutie = new Solution(solutieDeStart.getNumberOfVariables(), solutieDeStart.getNumberOfObjectives(), solutieDeStart.getNumberOfConstraints());
        for (int i = 0; i < numberOfVariables; i++) {
            RealVariable variable = (RealVariable) solutieDeStart.getVariable(i);
            double lowerBound = variable.getLowerBound();
            double upperBound = variable.getUpperBound();
            double valoareCoordonata = variable.getValue() + deviatiileStandard[i] * rand.nextGaussian();
            // corectieLaDepasireDomeniu
            if (valoareCoordonata > upperBound) {
                valoareCoordonata = upperBound;
            } else if (valoareCoordonata < lowerBound) {
                valoareCoordonata = lowerBound;
            }
            solutie.setVariable(i, new RealVariable(valoareCoordonata, lowerBound, upperBound));
        }
        return solutie;
    }

    private void calculeazaDeviatiileStandard(int index) {
        double suma;
        double leadersSize = leaders.size();
        for (int i = 0; i < numberOfVariables; i++) {
            suma = 0;
            for (int e = 0; e < leadersSize; e++) {
                double particle = EncodingUtils.getReal(particles[index].getVariable(1));
                double leader = EncodingUtils.getReal(leaders.get(e).getVariable(i));
                suma += abs(particle - leader);
                deviatiileStandard[i] = vitezaDeConvergenta * suma / leadersSize;
            }
        }
    }

}
