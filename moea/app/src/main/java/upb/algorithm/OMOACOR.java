package upb.algorithm;

import org.moeaframework.algorithm.pso.OMOPSO;
import org.moeaframework.analysis.sensitivity.EpsilonHelper;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.util.TypedProperties;

import java.util.Random;

// version 4
public class OMOACOR extends OMOPSO {

    private static final Random rand = new Random();
    private int numberOfVariables;
    private int numarFurnici;
    private final DEVIATIONS_CALCULATION_METHOD tipulCalcululuiDeviatiilor;
    private double vitezaDeConvergenta;
    private double[] deviatiileStandard;
    private int swarmSize;

    public OMOACOR(Problem problem, TypedProperties properties) {
        this(
                problem,
                properties.getInt("swarmSize", 100),
                properties.getInt("leaderSize", 100),
                new double[]{EpsilonHelper.getEpsilon(problem)},
                properties.getDouble("mutationProbability", 1.0 / problem.getNumberOfObjectives()),
                properties.getDouble("mutationPerturbation", 0.01),
                properties.getInt("maxIterations", 10000 / 100),
                properties.getDouble("convergeanceSpeed", 0.01),
                properties.getInt("numarFurnici", 100),
                DEVIATIONS_CALCULATION_METHOD.values()[properties.getInt("deviationsCalculationMethod", 7) - 1]
        );
    }

    public OMOACOR(
            Problem problem,
            int swarmSize,
            int leaderSize,
            double[] epsilons,
            double mutationProbability,
            double mutationPerturbation,
            int maxIterations,
            double vitezaDeConvergenta,
            int numarFurnici,
            DEVIATIONS_CALCULATION_METHOD tipulCalcululuiDeviatiilor
    ) {
        super(
                problem,
                swarmSize,
                leaderSize,
                epsilons,
                mutationProbability,
                mutationPerturbation,
                maxIterations
        );
        numberOfVariables = problem.getNumberOfVariables();
        this.vitezaDeConvergenta = vitezaDeConvergenta;
        deviatiileStandard = new double[swarmSize];
        this.numarFurnici = numarFurnici;
        this.tipulCalcululuiDeviatiilor = tipulCalcululuiDeviatiilor;
        this.swarmSize = swarmSize;
    }

    @Override
    protected void iterate() {

        try {
            constructieSolutii();
        } catch (Exception e) {
            return;
        }

        mutate();

        evaluateAll(particles);

        updateLocalBest();

        leaders.addAll(particles);
        leaders.update();

        if (archive != null) {
            archive.addAll(particles);
        }
    }

    private void constructieSolutii() throws Exception {
        for (int i = 0; i < particles.length; i++) {
            calculeazaDeviatiileStandard(i);
            particles[i] = construiesteSolutieNoua(i);
        }
    }

    private Solution construiesteSolutieNoua(int indexArhiva4Constructie) {
        Solution solutieDeStart = switch (tipulCalcululuiDeviatiilor) {
            case V1, V2, V3 -> particles[indexArhiva4Constructie];
            case V4, V5, V6 -> localBestParticles[indexArhiva4Constructie];
            default -> leaders.get(rand.nextInt(leaders.size()));
        };
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

    private void calculeazaDeviatiileStandard(int index) throws Exception {
        double suma;
        Solution leader;
        for (int i = 0; i < numberOfVariables; i++) {
            suma = 0;
            switch (tipulCalcululuiDeviatiilor) {
                case V1 -> {
                    for (int e = 0; e < swarmSize; e++)
                        suma = suma + Math.abs(((RealVariable) localBestParticles[e].getVariable(i)).getValue()
                                - ((RealVariable) particles[index].getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma / (swarmSize - 1);
                }
                case V2 -> {
                    leader = leaders.get(rand.nextInt(leaders.size()));
                    suma = Math.abs(((RealVariable) particles[index].getVariable(i)).getValue()
                            - ((RealVariable) leader.getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma;
                }
                case V3 -> {
                    for (int e = 0; e < leaders.size(); e++)
                        suma = suma + Math.abs(((RealVariable) particles[index].getVariable(i)).getValue()
                                - ((RealVariable) leaders.get(e).getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma / leaders.size();
                }
                case V4 -> {
                    for (int e = 0; e < swarmSize; e++)
                        suma = suma + Math.abs(((RealVariable) localBestParticles[e].getVariable(i)).getValue()
                                - ((RealVariable) localBestParticles[index].getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma / (swarmSize - 1);
                }
                case V5 -> {
                    leader = leaders.get(rand.nextInt(leaders.size()));
                    suma = Math.abs(((RealVariable) localBestParticles[index].getVariable(i)).getValue()
                            - ((RealVariable) leader.getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma;
                }
                case V6, V7 -> {
                    for (int e = 0; e < leaders.size(); e++)
                        suma = suma + Math.abs(((RealVariable) localBestParticles[index].getVariable(i)).getValue()
                                - ((RealVariable) leaders.get(e).getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma / leaders.size();
                }
            }
        }
    }

    public enum DEVIATIONS_CALCULATION_METHOD {
        V1, V2, V3, V4, V5, V6, V7
    }

}
