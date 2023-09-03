package upb;

import org.moeaframework.core.spi.AlgorithmFactory;
import org.moeaframework.core.spi.ProblemFactory;
import upb.algorithm.ACORProvider;
import upb.problem.LoneyMOProvider;

public class Main {
    public static void main(String[] args) {
        AlgorithmFactory.getInstance().addProvider(new ACORProvider());
        ProblemFactory.getInstance().addProvider(new LoneyMOProvider());
    }
}
