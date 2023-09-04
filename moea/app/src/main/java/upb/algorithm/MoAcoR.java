package upb.algorithm;

import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.algorithm.ReferencePointNondominatedSortingPopulation;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.AggregateConstraintComparator;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.DominanceComparator;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.core.operator.TournamentSelection;
import org.moeaframework.core.spi.OperatorFactory;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.util.TypedProperties;
import org.moeaframework.util.weights.NormalBoundaryDivisions;

import java.util.Random;

// version 1
public class MoAcoR extends AbstractEvolutionaryAlgorithm implements EpsilonBoxEvolutionaryAlgorithm {

    private final TypedProperties properties = new TypedProperties();

    public MoAcoR(Problem problem, TypedProperties typedProperties) {
        super(
                problem,
                100,
                new ReferencePointNondominatedSortingPopulation(
                        problem.getNumberOfObjectives(),
                        new NormalBoundaryDivisions(1, 1)
                ),
                null,
                new RandomInitialization(problem),
                OperatorFactory.getInstance().getVariation(problem)
        );

        properties.setDouble("localitateaProcesuluiDeCautare", 0.1);
        properties.setDouble("vitezaDeConvergenta", 0.85);
        properties.setInt("numarFurnici", 100);
        properties.setBoolean("selectieTournamet", false);

    }

    @Override
    public void iterate() {
        NondominatedSortingPopulation population = getPopulation();
        EpsilonBoxDominanceArchive archive = getArchive();
        Problem problem = getProblem();

        MoAcoRUtils moAcoRUtils = new MoAcoRUtils(problem, properties, population, this.getVariation());
        Population offspring = moAcoRUtils.constructieSolutii_Evaluare_CalculareFurnicaRBestGlobala();

        evaluateAll(offspring);

        if (archive != null) {
            archive.addAll(offspring);
        }

        population.addAll(offspring);
        population.truncate(offspring.size());
    }

    @Override
    public EpsilonBoxDominanceArchive getArchive() {
        return (EpsilonBoxDominanceArchive) super.getArchive();
    }

    @Override
    public NondominatedSortingPopulation getPopulation() {
        return (NondominatedSortingPopulation) super.getPopulation();
    }

    private static class MoAcoRUtils {

        private final Population offspring;
        private final int populationSize;
        private final int numberOfVariables;
        private final int numarFurnici;
        private final TypedProperties properties;
        private final Random rand;
        private final Population population;
        private final MoAcoRTournamentSelection selection;
        private final Variation variation;

        private final double[] deviatiileStandard;
        private final double[] ponderiW;
        private final double[] probabilitatiDeSelectie;
        private final double[] probabilitatiCumulateDeSelectie;

        public MoAcoRUtils(Problem problem, TypedProperties properties, Population population, Variation variation) {
            this.numberOfVariables = problem.getNumberOfVariables();
            this.properties = properties;
            this.population = population;
            this.offspring = new Population();
            this.populationSize = population.size();
            this.numarFurnici = properties.getInt("numarFurnici", 10);
            this.variation = variation;
            rand = new Random();
            selection = new MoAcoRTournamentSelection();
            deviatiileStandard = new double[numberOfVariables];
            ponderiW = new double[populationSize];
            probabilitatiDeSelectie = new double[populationSize];
            probabilitatiCumulateDeSelectie = new double[populationSize];
            calculeazaPonderiSiProbabilitatiDeSelectie();
        }

        public Population constructieSolutii_Evaluare_CalculareFurnicaRBestGlobala() {
            constructieSolutii();
            return offspring;
        }

        private void constructieSolutii() {
            while (offspring.size() < numarFurnici) {
                int indexArhiva4Constructie;
                if (properties.getBoolean("selectieTournamet", false)) {
                    indexArhiva4Constructie = selection.selectIndex(population);
                } else {
                    indexArhiva4Constructie = selectieIndexArhiva4Constructie();
                }

                calculeazadeviatiileStandard(indexArhiva4Constructie);
                Solution sol1 = construiesteSolutieNoua(indexArhiva4Constructie);
                if (properties.getBoolean("selectieTournamet", false)) {
                    indexArhiva4Constructie = selection.selectIndex(population);
                } else {
                    indexArhiva4Constructie = selectieIndexArhiva4Constructie();
                }

                calculeazadeviatiileStandard(indexArhiva4Constructie);
                Solution sol2 = construiesteSolutieNoua(indexArhiva4Constructie);
                offspring.addAll(variation.evolve(new Solution[]{sol1, sol2}));
            }
        }

        private void calculeazaPonderiSiProbabilitatiDeSelectie() {
            double constanta = properties.getDouble("localitateaProcesuluiDeCautare", 1e-3) * populationSize * Math.sqrt(2 * Math.PI);// q k * radical (2 * PI);
            for (int i = 0; i < populationSize; i++) {
                ponderiW[i] = Math.exp((-1) * Math.pow(i, 2) / (2 * Math.pow(constanta, 2) * Math.pow(populationSize, 2))) / constanta;
            }
            // calculeaza probabilitati de selectie (inclusiv cumulate)
            double sumaPonderi = 0;
            for (int i = 0; i < populationSize; i++) {
                sumaPonderi += ponderiW[i];
            }
            for (int i = 0; i < populationSize; i++) {
                probabilitatiDeSelectie[i] = ponderiW[i] / sumaPonderi;
            }

            probabilitatiCumulateDeSelectie[0] = probabilitatiDeSelectie[0];
            probabilitatiCumulateDeSelectie[populationSize - 1] = 1;
            for (int i = 1; i < populationSize - 1; i++) {
                probabilitatiCumulateDeSelectie[i] = probabilitatiCumulateDeSelectie[i - 1] + probabilitatiDeSelectie[i];
            }
        }

        private int selectieIndexArhiva4Constructie() {
            // selectie ruleta
            double numarAleator = rand.nextDouble();
            for (int i = 0; i < populationSize; i++) {
                if (numarAleator < probabilitatiCumulateDeSelectie[i]) {
                    return i;
                }
            }
            // shoould NEVER be here
            return -1;
        }

        private void calculeazadeviatiileStandard(int indexArhiva4Constructie) {
            double suma;
            for (int i = 0; i < numberOfVariables; i++) {
                suma = 0;
                for (int e = 0; e < populationSize; e++)
                    suma = suma + Math.abs(((RealVariable) population.get(e).getVariable(i)).getValue() - ((RealVariable) population.get(indexArhiva4Constructie).getVariable(i)).getValue());
                deviatiileStandard[i] = properties.getDouble("vitezaDeConvergenta", 0.85) * suma / (populationSize - 1);
            }
        }

        private Solution construiesteSolutieNoua(int indexArhiva4Constructie) {
            Solution solutieDeStart = population.get(indexArhiva4Constructie);
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

        private static class MoAcoRTournamentSelection extends TournamentSelection {

            private final int size;

            public MoAcoRTournamentSelection() {
                super(2, new ChainedComparator(new AggregateConstraintComparator(), (solution1, solution2) -> PRNG.nextBoolean() ? -1 : 1));
                this.size = 2;
            }

            public int selectIndex(Population population) {
                DominanceComparator comparator = getComparator();
                int winnerIndex = PRNG.nextInt(population.size());
                Solution winner = population.get(winnerIndex);

                for (int i = 1; i < size; i++) {
                    int candidateIndex = PRNG.nextInt(population.size());
                    Solution candidate = population.get(candidateIndex);

                    int flag = comparator.compare(winner, candidate);

                    if (flag > 0) {
                        winner = candidate;
                        winnerIndex = candidateIndex;
                    }
                }

                return winnerIndex;
            }

        }
    }

}