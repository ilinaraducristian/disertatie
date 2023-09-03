package upb.problem;

import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Problem;
import org.moeaframework.core.spi.ProblemProvider;

import java.io.IOException;

import static org.moeaframework.core.PopulationIO.readReferenceSet;

public class LoneyMOProvider extends ProblemProvider {
    @Override
    public Problem getProblem(String name) {
        if (!name.equals("LoneyMO")) {
            return null;
        }
        return new LoneyMO();
    }

    @Override
    public NondominatedPopulation getReferenceSet(String name) {
        if (!name.equals("LoneyMO")) {
            return null;
        }
        try {
            return readReferenceSet("LoneyMO");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
