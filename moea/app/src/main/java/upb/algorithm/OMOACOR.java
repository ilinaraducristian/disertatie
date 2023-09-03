package upb.algorithm;

import org.moeaframework.algorithm.pso.OMOPSO;
import org.moeaframework.analysis.sensitivity.EpsilonHelper;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.RealVariable;

import java.util.Random;

// version 7
public class OMOACOR extends OMOPSO {

    private static final Random rand = new Random();
    int numberOfVariables;
    int numarFurnici;
    int tipulCalcululuiDeviatiilor;
    double vitezaDeConvergenta;
    double[] deviatiileStandard;

    int swarmSize;

    public OMOACOR(Problem problem) {
        this(problem, 100, 100,
                new double[]{EpsilonHelper.getEpsilon(problem)}, 1.0 / problem.getNumberOfObjectives(), 0.01, 10000 / 100,
                0.01f, 100, 7);
    }

    public OMOACOR(Problem problem, double[] epsilon, double mutationProbability, double mutationPerturbation, double convergenceSpeed) {
        this(problem, 100, 100, epsilon, mutationProbability, mutationPerturbation, 10000 / 100,
                convergenceSpeed, 100, 7);
    }

    public OMOACOR(Problem problem,
                   int swarmSize,
                   int leaderSize,
                   double[] epsilons,
                   double mutationProbability,
                   double mutationPerturbation,
                   int maxIterations,
                   double vitezaDeConvergenta,
                   int numarFurnici,
                   int tipulCalcululuiDeviatiilor) {
        super(problem, swarmSize, leaderSize, epsilons, mutationProbability, mutationPerturbation, maxIterations);
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
            case 1, 2, 3 -> particles[indexArhiva4Constructie];
            case 4, 5, 6 -> localBestParticles[indexArhiva4Constructie];
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
        Solution leader = null;
        for (int i = 0; i < numberOfVariables; i++) {
            suma = 0;
            switch (tipulCalcululuiDeviatiilor) {
                case 1:
                    for (int e = 0; e < swarmSize; e++)
                        suma = suma + Math.abs(((RealVariable) localBestParticles[e].getVariable(i)).getValue()
                                - ((RealVariable) particles[index].getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma / (swarmSize - 1);
                    break;
                case 2:
                    leader = leaders.get(rand.nextInt(leaders.size()));
                    suma = Math.abs(((RealVariable) particles[index].getVariable(i)).getValue()
                            - ((RealVariable) leader.getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma;
                    break;
                case 3:
                    for (int e = 0; e < leaders.size(); e++)
                        suma = suma + Math.abs(((RealVariable) particles[index].getVariable(i)).getValue()
                                - ((RealVariable) leaders.get(e).getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma / leaders.size();
                    break;
                case 4:
                    for (int e = 0; e < swarmSize; e++)
                        suma = suma + Math.abs(((RealVariable) localBestParticles[e].getVariable(i)).getValue()
                                - ((RealVariable) localBestParticles[index].getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma / (swarmSize - 1);
                    break;
                case 5:
                    leader = leaders.get(rand.nextInt(leaders.size()));
                    suma = Math.abs(((RealVariable) localBestParticles[index].getVariable(i)).getValue()
                            - ((RealVariable) leader.getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma;
                    break;
                case 6:
                case 7:
                    for (int e = 0; e < leaders.size(); e++)
                        suma = suma + Math.abs(((RealVariable) localBestParticles[index].getVariable(i)).getValue()
                                - ((RealVariable) leaders.get(e).getVariable(i)).getValue());
                    deviatiileStandard[i] = vitezaDeConvergenta * suma / leaders.size();
                    break;
                default:
                    throw new Exception("tipulCalcululuiDeviatiilor cannot be 0");
            }
        }
    }

}
